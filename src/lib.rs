pub mod bai_parser;
pub mod bam_parser;
pub mod timer;
use bstr::{BString, ByteSlice};
use noodles::bam::io::reader::header;
use noodles::bam::record;
use noodles::bam::{self as noodles_bam, record::Record};
use noodles::bgzf::{VirtualPosition, io as bgzf_io};
use noodles::sam::{
    self as noodles_sam, Header,
    alignment::{RecordBuf, io::Write},
};
use rayon::prelude::*;
use std::fs::File;
use std::{
    collections::HashMap,
    io,
    num::NonZero,
    sync::{Arc, Mutex},
};

//TODO: move count records to here
//TODO: implement BAM storm func(F)
//TODO: How to write results to BAM file?

// Get BAM header
pub fn get_bam_header(bam_path: &str) -> io::Result<noodles_sam::Header> {
    let bam_file = File::open(bam_path)?;
    let mut reader = noodles_bam::io::Reader::new(bam_file);
    let header = reader.read_header()?;
    Ok(header)
}

// Count all records from a standard Noodles BAM reader
pub fn count_from_standard_bam_reader(bam_path: &str, thread_num: u64) -> io::Result<u64> {
    // Set Single/Multiple threads Reader
    let worker_count = NonZero::new(thread_num as usize).unwrap_or(NonZero::<usize>::MIN);
    // Set Multithread Reader
    let bam_file = File::open(bam_path)?;
    let mut reader = noodles_bam::io::Reader::from(
        bgzf_io::MultithreadedReader::with_worker_count(worker_count, bam_file),
    );
    let _header = reader.read_header()?;
    let mut n_records: u64 = 0;
    for result in reader.records() {
        let _record = result?;
        n_records += 1;
        // Process the record as needed
        // For demonstration, we will just print the read name
    }
    println!(
        "===== Total number of records from standard reader: {}.",
        n_records
    );
    Ok(n_records)
}

// Count all records from a BAM file
pub fn count_all_records(
    bam_path: &str,
    intervals: &[(VirtualPosition, VirtualPosition)],
    thread_num: u64,
) -> io::Result<u64> {
    let all_interval_readers = bam_parser::get_entire_bam_reader(bam_path, intervals, thread_num)?;
    let total_records: u64 = all_interval_readers
        .into_par_iter()
        .map(|mut reader| {
            let mut local_count: u64 = 0;
            for result in reader.records() {
                let _record = result?;
                local_count += 1;
            }
            Ok(local_count)
        })
        .sum::<io::Result<u64>>()?;

    println!(
        "===== Total number of records from interval reader: {}.",
        total_records
    );
    Ok(total_records)
}

// Process records from BAM file in parallel through intervals
pub fn process_records_through_intervals<F>(
    bam_path: &str,
    intervals: &[(VirtualPosition, VirtualPosition)],
    thread_num: u64,
    output_writer: Arc<Mutex<noodles_bam::io::Writer<bgzf_io::Writer<File>>>>,
    process_record: F,
) -> io::Result<()>
where
    F: Fn(
            &noodles_bam::record::Record,
            &Arc<Mutex<noodles_bam::io::Writer<bgzf_io::Writer<File>>>>,
        ) -> Result<(), Box<dyn std::error::Error>>
        + Send
        + Sync,
{
    let all_interval_readers = bam_parser::get_entire_bam_reader(bam_path, intervals, thread_num)?;

    all_interval_readers
        .into_par_iter()
        .try_for_each(|mut reader| -> io::Result<()> {
            for result in reader.records() {
                let record = result?;
                process_record(&record, &output_writer).unwrap();
            }
            Ok(())
        })?;

    Ok(())
}

fn replace_phi(
    original_string: &str,
    phi_pattern: &str,
    new_string: &str,
) -> Result<String, std::io::Error> {
    let updated_string = original_string.replace(phi_pattern, new_string);
    Ok(updated_string)
}

// Convert a generic SAM Record/RecordBuf into a BAM record
pub fn try_into_bam_record<R>(
    header: &noodles_sam::Header,
    record: &R,
) -> io::Result<noodles_bam::Record>
where
    R: noodles_sam::alignment::Record,
{
    let mut writer = noodles_bam::io::Writer::from(Vec::new());
    writer.write_alignment_record(header, record)?;

    let src = writer.into_inner();
    let mut reader = noodles_bam::io::Reader::from(&src[..]);
    let mut record = noodles_bam::Record::default();
    reader.read_record(&mut record)?;

    Ok(record)
}

pub fn update_bam_record(
    header: &Header,
    record: &Record,
    updated_fields: &HashMap<&str, &str>,
    phi_pattern: &str,
    output_writer: &Arc<Mutex<noodles_bam::io::Writer<bgzf_io::Writer<File>>>>,
) -> io::Result<()> {
    let mut updated_record = RecordBuf::try_from_alignment_record(header, record).unwrap();

    // Update the record fields based on the provided HashMap
    for (&field_name, &field_value) in updated_fields {
        match field_name {
            "name" => {
                // Update the name field with the new value
                let name_pointer = updated_record.name_mut();
                let original_string = record.name().unwrap();
                // convert
                let updated_string = original_string
                    .to_str()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                    .and_then(|s| replace_phi(s, phi_pattern, field_value))?;
                *name_pointer = Some(BString::from(updated_string));
            }
            _ => {
                eprintln!("Field {} not found", field_name);
            }
        }
    }
    let updated_record = try_into_bam_record(header, &updated_record)?;
    output_writer
        .lock()
        .unwrap()
        .write_record(&header, &updated_record)?;
    Ok(())
}

// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bai_parser::{get_linear_indexes, get_linear_intervals};
    use crate::bam_parser::count_from_standard_bam_reader;
    use crate::timer::timeit;

    #[test]
    fn test_count_through_intervals() -> Result<(), Box<dyn std::error::Error>> {
        let test_bam = "tests/chr_all.bam";
        let bai_path = "tests/chr_all.bam.bai";
        println!(
            "===== Testing read through intervals for BAM: {:?}.",
            test_bam
        );
        let linear_indexes_all = get_linear_indexes(bai_path)?;
        let intervals = get_linear_intervals(&linear_indexes_all)?;
        let standard_reader_total_records: u64 = count_from_standard_bam_reader(test_bam, 4)?;
        let all_intervals_total_records = timeit(|| count_all_records(test_bam, &intervals, 4))?;
        assert_eq!(all_intervals_total_records, standard_reader_total_records);
        Ok(())
    }
}
