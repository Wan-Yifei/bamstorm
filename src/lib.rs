pub mod bai_parser;
pub mod bam_parser;
pub mod timer;
#[cfg(all(feature = "python", unix))]
#[allow(unsafe_op_in_unsafe_fn)]
pub mod python;
use bstr::{BString, ByteSlice};
use noodles::bam::{self as noodles_bam, record::Record};
use noodles::bgzf::{VirtualPosition, io as bgzf_io};
use noodles::sam::{
    self as noodles_sam, Header,
    alignment::{RecordBuf, io::Write as _},
};
use rayon::prelude::*;
use std::fs::File;
use std::{
    collections::HashMap,
    io::{self, Write},
    num::NonZero,
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
) -> io::Result<u64> {
    let all_intervals = bam_parser::get_entire_bam_intervals(bam_path, intervals)?;
    let all_intervals = bam_parser::merge_intervals(&all_intervals, rayon::current_num_threads());
    let total_records: u64 = all_intervals
        .into_par_iter()
        .map(|(start, end)| -> io::Result<u64> {
            let mut reader = bam_parser::read_bam_by_interval(bam_path, start, end)?;
            let mut count = 0u64;
            for result in reader.records() {
                result?;
                count += 1;
            }
            Ok(count)
        })
        .sum::<io::Result<u64>>()?;

    println!(
        "===== Total number of records from interval reader: {}.",
        total_records
    );
    Ok(total_records)
}

// BGZF EOF block (28 bytes) — stripped from intermediate chunks, appended once at end
const BGZF_EOF: &[u8; 28] = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\
    \x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";

// Process records from BAM file in parallel through intervals.
// Each rayon worker writes to a local in-memory BGZF buffer; rayon collect()
// preserves interval order, so the concatenated output is coordinate-sorted.
pub fn process_records_through_intervals<F>(
    bam_path: &str,
    intervals: &[(VirtualPosition, VirtualPosition)],
    header: &noodles_sam::Header,
    output_writer: noodles_bam::io::Writer<bgzf_io::Writer<File>>,
    process_record: F,
) -> io::Result<()>
where
    F: Fn(&noodles_bam::record::Record) -> io::Result<Option<RecordBuf>> + Send + Sync,
{
    let all_intervals = bam_parser::get_entire_bam_intervals(bam_path, intervals)?;

    // IO and processing are fused per-interval: each rayon worker reads its slice
    // from disk, processes records, and compresses output independently.
    // Concurrency is bounded by rayon's thread pool (≈ CPU cores), so at most
    // thread_count interval buffers are live simultaneously — O(cores × interval_size)
    // peak memory instead of O(entire file). The tail interval is also included here.
    let chunks: Vec<Vec<u8>> = all_intervals
        .into_par_iter()
        .map(|(start, end)| -> io::Result<Vec<u8>> {
            let mut reader = bam_parser::read_bam_by_interval(bam_path, start, end)?;
            let mut local_writer =
                noodles_bam::io::Writer::from(bgzf_io::Writer::new(Vec::new()));
            for result in reader.records() {
                let record = result?;
                if let Some(processed) = process_record(&record)? {
                    local_writer.write_alignment_record(header, &processed)?;
                }
            }
            local_writer.into_inner().finish()
        })
        .collect::<io::Result<Vec<_>>>()?;

    // Consume the writer chain to reach the raw File.
    // flush() writes any pending uncompressed bytes as a BGZF block (no EOF).
    // into_inner() on the bgzf Writer takes the File; Drop sees inner=None and
    // skips try_finish(), so no spurious EOF block is written between header and records.
    let mut bgzf_writer = output_writer.into_inner();
    bgzf_writer.flush()?;
    let mut file = bgzf_writer.into_inner();

    // Strip EOF from every chunk and write data blocks in order,
    // then append a single EOF block to close the BGZF stream.
    for chunk in &chunks {
        let data_end = chunk.len().saturating_sub(BGZF_EOF.len());
        file.write_all(&chunk[..data_end])?;
    }
    file.write_all(BGZF_EOF)?;

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

pub fn update_bam_record(
    header: &Header,
    record: &Record,
    updated_fields: &HashMap<&str, &str>,
    phi_pattern: &str,
) -> io::Result<RecordBuf> {
    let mut updated_record = RecordBuf::try_from_alignment_record(header, record)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    for (&field_name, &field_value) in updated_fields {
        match field_name {
            "name" => {
                let original_string = record
                    .name()
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "record missing name"))?;
                let updated_string = original_string
                    .to_str()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                    .and_then(|s| replace_phi(s, phi_pattern, field_value))?;
                *updated_record.name_mut() = Some(BString::from(updated_string));
            }
            _ => {
                eprintln!("Field {} not found", field_name);
            }
        }
    }
    Ok(updated_record)
}

// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bai_parser::{get_linear_indexes, get_linear_intervals};
    use crate::timer::timeit;
    use noodles::bam as noodles_bam;
    use std::collections::HashMap;

    const TEST_BAM: &str = "tests/mt.sorted.bam";
    const TEST_BAI: &str = "tests/mt.sorted.bam.bai";

    #[test]
    fn test_count_through_intervals() -> Result<(), Box<dyn std::error::Error>> {
        println!(
            "===== Testing read through intervals for BAM: {:?}.",
            TEST_BAM
        );
        let linear_indexes_all = get_linear_indexes(TEST_BAI)?;
        let intervals = get_linear_intervals(&linear_indexes_all)?;
        let standard_reader_total_records: u64 =
            crate::count_from_standard_bam_reader(TEST_BAM, 4)?;
        let all_intervals_total_records =
            timeit(|| count_all_records(TEST_BAM, &intervals))?;
        assert_eq!(all_intervals_total_records, standard_reader_total_records);
        Ok(())
    }

    // Verify that update_bam_record correctly replaces a substring in read names
    // and that the output BAM can be read back with the expected modification.
    #[test]
    fn test_update_bam_record_name_substitution() -> Result<(), Box<dyn std::error::Error>> {
        let output_path = "tests/test_output_name_sub.bam";
        let header = get_bam_header(TEST_BAM)?;

        let mut writer = noodles_bam::io::Writer::new(std::fs::File::create(output_path)?);
        writer.write_header(&header)?;

        let mut updated_fields = HashMap::new();
        updated_fields.insert("name", "REPLACED");

        let linear_indexes_all = get_linear_indexes(TEST_BAI)?;
        let intervals = get_linear_intervals(&linear_indexes_all)?;

        let phi_pattern = "-";
        process_records_through_intervals(
            TEST_BAM,
            &intervals,
            &header,
            writer,
            |record| {
                let updated =
                    update_bam_record(&header, record, &updated_fields, phi_pattern)?;
                Ok(Some(updated))
            },
        )?;

        // Read back and verify: every read name should contain "REPLACED" (replacing "-")
        let out_file = std::fs::File::open(output_path)?;
        let mut reader = noodles_bam::io::Reader::new(out_file);
        let _hdr = reader.read_header()?;
        let mut checked = 0u64;
        for result in reader.records() {
            let record = result?;
            if let Some(name) = record.name() {
                let name_str = name.to_str().unwrap_or("");
                assert!(
                    name_str.contains("REPLACED"),
                    "Expected 'REPLACED' in name, got: {name_str}"
                );
                checked += 1;
                if checked >= 20 {
                    break;
                }
            }
        }
        assert!(checked > 0, "No records found in output BAM");
        std::fs::remove_file(output_path)?;
        Ok(())
    }

    // Verify process_records_through_intervals output record count equals standard reader.
    // This exercises the lazy IO path (Step 5) and proves the tail interval is included
    // in the parallel pass (Step 6) — a missing tail interval would produce a lower count.
    #[test]
    fn test_process_intervals_record_count() -> Result<(), Box<dyn std::error::Error>> {
        let output_path = "tests/test_output_count.bam";
        let header = get_bam_header(TEST_BAM)?;
        let standard_count = crate::count_from_standard_bam_reader(TEST_BAM, 1)?;

        let mut writer = noodles_bam::io::Writer::new(std::fs::File::create(output_path)?);
        writer.write_header(&header)?;

        let linear_indexes_all = get_linear_indexes(TEST_BAI)?;
        let intervals = get_linear_intervals(&linear_indexes_all)?;

        process_records_through_intervals(
            TEST_BAM,
            &intervals,
            &header,
            writer,
            |record| Ok(Some(RecordBuf::try_from_alignment_record(&header, record)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?)),
        )?;

        // Count records in output and compare with standard reader
        let out_file = std::fs::File::open(output_path)?;
        let mut reader = noodles_bam::io::Reader::new(out_file);
        let _hdr = reader.read_header()?;
        let output_count = reader.records().count() as u64;

        assert_eq!(
            output_count, standard_count,
            "process_records_through_intervals output {output_count} != standard reader {standard_count}"
        );
        std::fs::remove_file(output_path)?;
        Ok(())
    }
}
