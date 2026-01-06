use noodles::bam;
use noodles::bgzf::{VirtualPosition, io as bgzf_io};
use rayon::prelude::*;
use std::{
    fs::File,
    io::{self, Cursor, Seek as f_seek, SeekFrom, prelude::*},
    num::NonZero,
};

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
    bam_file.read_exact(&mut buffer).unwrap_or_else(|e| {
        panic!(
            "read_exact failed: start_voffset={:?}, end_voffset={:?}, \
         compressed_range=[{}..{}], buffer_size={} bytes, error={}",
            start_voffset,
            end_voffset,
            start_voffset.compressed(),
            end_voffset.compressed(),
            buffer_size,
            e
        )
    });
    // Suppose `buffer` is Vec<u8> read from the BAM
    let bam_interval_cursor = Cursor::new(buffer);
    let decoder =
        bgzf_io::MultithreadedReader::with_worker_count(worker_count, bam_interval_cursor);
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
            read_bam_by_interval(bam_path, start_voffset, end_voffset, thread_num)
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
        bam_file.read_to_end(&mut end_buffer).unwrap_or_else(|e| {
            panic!(
                "read_to_end failed: final voffset {:?}, final_interval_coffset={}, \
             error={}",
                last_end_voffset, final_interval_coffset, e
            )
        });
        let final_interval_cursor = Cursor::new(end_buffer);
        let worker_count = NonZero::new(thread_num as usize).unwrap_or(NonZero::<usize>::MIN);
        let decoder =
            bgzf_io::MultithreadedReader::with_worker_count(worker_count, final_interval_cursor);
        let reader_eof = bam::io::Reader::from(decoder);
        all_interval_readers.push(reader_eof);
    }

    Ok(all_interval_readers)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::bai_parser::{get_linear_indexes, get_linear_intervals};
    use crate::timer::timeit;

    #[test]
    #[ignore]
    fn test_bam_read_by_interval() -> Result<(), Box<dyn std::error::Error>> {
        let test_bam = "tests/chr_all.bam";
        let bai_path = "tests/chr_all.bam.bai";
        // let bai_path = "tests/full.bam.bai";
        // let test_bam = "tests/full.bam";
        let linear_indexes_all = get_linear_indexes(bai_path)?;
        let intervals = get_linear_intervals(&linear_indexes_all)?;
        timeit(|| read_bam_by_interval(test_bam, intervals[0].0, intervals[0].1, 4))?;
        // println!("{:?}", test_bam);
        Ok(())
    }
}
