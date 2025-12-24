use noodles::bam;
use noodles::bgzf::{VirtualPosition, io as bgzf_io};
use std::{
    fs::File,
    io::{self, prelude::*, Seek as f_seek, SeekFrom, Cursor},
    num::NonZero,
};
use rayon::prelude::*;


pub fn standard_bam_read(bam_path: &str, thread_num: u64) -> io::Result<u64> {
    // Set Single/Multiple threads Reader
    let worker_count = NonZero::new(thread_num as usize).unwrap_or(NonZero::<usize>::MIN);
    // Set Multithread Reader
    let bam_file = File::open(bam_path)?;
    let mut reader = bam::io::Reader::from(
        bgzf_io::MultithreadedReader::with_worker_count(worker_count, bam_file)
    );
    let _header = reader.read_header()?;
    let mut n_records: u64 = 0;
    for result in reader.records() {
        let _record = result?;
        n_records += 1;
        // Process the record as needed
        // For demonstration, we will just print the read name
    }   
    println!("===== Total number of records from standard reader: {}.", n_records);
    Ok(n_records)
}


// TODO: implement reading BAM by intervals
// TODO: 1. Seek to the start virtual position
// TODO: 2. Read until reaching the end virtual position
// TODO: 3. multithreading support for reading BAM by intervals
// TODO: 4. Enable multithreading decompression of BGZF blocks
// TODO: 5. Need a method to append the vector of linear indexes to point to the end of the BAM file 
pub fn read_bam_by_interval(
    bam_path: &str,
    start_voffset: VirtualPosition,
    end_voffset: VirtualPosition,
    thread_num: u64,
) -> io::Result<bam::io::Reader<bgzf_io::MultithreadedReader<Cursor<Vec<u8>>>>> {
    // Set threads
    // let worker_count = thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN);
    let worker_count = NonZero::new(thread_num as usize).unwrap_or(NonZero::<usize>::MIN);
    // Set Multithread Reader
    let mut bam_file = File::open(bam_path)?;
    // Seek to the start coffset of the interval
    bam_file.seek(SeekFrom::Start(start_voffset.compressed()))?;
    // println!("coffset {:?} of voffset: {:?}.", start_voffset.compressed(), start_voffset);
    let buffer_size = end_voffset.compressed() - start_voffset.compressed(); // size to read complete BGZF blocks 
    let mut buffer = vec![0; buffer_size as usize];
    // println!("Reading interval buffer size: {:?} bytes.", buffer_size);
    bam_file.read_exact(&mut buffer)?;
    // Suppose `buffer` is Vec<u8> read from the BAM
    let bam_interval_cursor = Cursor::new(buffer);
    let decoder = bgzf_io::MultithreadedReader::with_worker_count(worker_count, bam_interval_cursor);
    let reader = bam::io::Reader::from(decoder);
    Ok(reader)
}

pub fn read_through_intervals(
    bam_path: &str,
    intervals: &[(VirtualPosition, VirtualPosition)],
    thread_num: u64,
) -> io::Result<Vec<bam::io::Reader<bgzf_io::MultithreadedReader<Cursor<Vec<u8>>>>>> {

    // Run intervals in parallel
    let all_interval_readers = intervals
        .par_iter()
        .map(|&(start_voffset, end_voffset)| {
            read_bam_by_interval(
                bam_path,
                start_voffset,
                end_voffset,
                thread_num,
            )
        })
        .collect::<Result<Vec<_>, io::Error>>()?;

        Ok(all_interval_readers)
}

pub fn get_entire_bam_reader(
    bam_path: &str,
    intervals: &[(VirtualPosition, VirtualPosition)],
    thread_num: u64,
) -> io::Result<Vec<bam::io::Reader<bgzf_io::MultithreadedReader<Cursor<Vec<u8>>>>>> {

    let mut all_interval_readers = read_through_intervals(bam_path, intervals, thread_num)?;

    // Add the final interval: [last end_coffset, EOF]
    if let Some((_, last_end_voffset)) = intervals.last() {
        let mut bam_file = File::open(bam_path)?;
        let final_interval_coffset = last_end_voffset.compressed();
        bam_file.seek(SeekFrom::Start(final_interval_coffset))?;
        let mut end_buffer: Vec<u8> = Vec::new();
        bam_file.read_to_end(&mut end_buffer)?;
        let final_interval_cursor = Cursor::new(end_buffer);
        let worker_count = NonZero::new(thread_num as usize).unwrap_or(NonZero::<usize>::MIN);
        let decoder = bgzf_io::MultithreadedReader::with_worker_count(worker_count, final_interval_cursor);
        let reader_eof = bam::io::Reader::from(decoder);
        all_interval_readers.push(reader_eof);
    }

    Ok(all_interval_readers)
}

pub fn count_all_records(
    bam_path: &str,
    intervals: &[(VirtualPosition, VirtualPosition)],
    thread_num: u64,
) -> io::Result<u64> {

    let mut total_records: u64 = 0;
    let all_interval_readers = get_entire_bam_reader(bam_path, intervals, thread_num)?;

    for mut reader in all_interval_readers {
        for result in reader.records() {
            let _record = result?;
            total_records += 1;
        }
    }

    println!("===== Total number of records from interval reader: {}.", total_records);
    Ok(total_records)
}


#[cfg(test)]
mod test {
    use crate::timer::timeit;
    use crate::bai_parser::{get_linear_indexes, get_linear_intervals};
    use super::*;

    #[test]
    fn test_standard_bam_read_with_timer() -> Result<(), Box<dyn std::error::Error>> {
        // Run bam_read and propagate any errors
        let test_bam = "/Users/yifeiwan/Projects/bamstorm/tests/chr1.bam";
        // let test_bam = "/Users/yifeiwan/Projects/bamstorm/tests/full.bam";
        timeit(|| standard_bam_read(test_bam, 4))?;
        Ok(())
    }

    #[test]
    #[ignore]
    fn test_bam_read_by_interval() -> Result<(), Box<dyn std::error::Error>> {
        let test_bam = "/Users/yifeiwan/Projects/bamstorm/tests/chr1.bam";
        let bai_path = "/Users/yifeiwan/Projects/bamstorm/tests/chr1.bam.bai";
        // let bai_path = "/Users/yifeiwan/Projects/bamstorm/full.bam.bai";
        // let test_bam = "/Users/yifeiwan/Projects/bamstorm/full.bam";
        let linear_indexes_all = get_linear_indexes(bai_path)?;
        let intervals = get_linear_intervals(&linear_indexes_all)?;
        timeit(|| read_bam_by_interval(test_bam, intervals[0].0, intervals[0].1, 4))?;
        // println!("{:?}", test_bam);
        Ok(())
    }

    #[test]
    fn test_read_through_intervals() -> Result<(), Box<dyn std::error::Error>> {
        // let test_bam = "/Users/yifeiwan/Projects/bamstorm/tests/full.bam";
        // let bai_path = "/Users/yifeiwan/Projects/bamstorm/tests/full.bam.bai";
        let test_bam = "/Users/yifeiwan/Projects/bamstorm/tests/chr1.bam";
        let bai_path = "/Users/yifeiwan/Projects/bamstorm/tests/chr1.bam.bai";
        let linear_indexes_all = get_linear_indexes(bai_path)?;
        let intervals = get_linear_intervals(&linear_indexes_all)?;
        let standard_reader_total_records:u64= standard_bam_read(test_bam, 4)?;
        let all_intervals_total_records = timeit(|| count_all_records(test_bam, &intervals, 4))?;
        assert_eq!(all_intervals_total_records, standard_reader_total_records);
        Ok(())
    }
 
}
