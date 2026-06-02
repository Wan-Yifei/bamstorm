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
) -> io::Result<bam::io::Reader<bgzf_io::MultithreadedReader<Cursor<Vec<u8>>>>> {
    let mut bam_file = File::open(bam_path)?;
    bam_file.seek(SeekFrom::Start(start_voffset.compressed()))?;
    let buffer_size = end_voffset.compressed() - start_voffset.compressed();
    let mut buffer = vec![0; buffer_size as usize];
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
    // Each interval reader uses a single decompression worker; outer rayon parallelism
    // across intervals is sufficient — adding more workers per reader over-subscribes the CPU.
    let decoder = bgzf_io::MultithreadedReader::with_worker_count(
        NonZero::<usize>::MIN,
        Cursor::new(buffer),
    );
    Ok(bam::io::Reader::from(decoder))
}

pub fn read_through_intervals(
    bam_path: &str,
    intervals: &[(VirtualPosition, VirtualPosition)],
) -> io::Result<Vec<bam::io::Reader<bgzf_io::MultithreadedReader<Cursor<Vec<u8>>>>>> {
    intervals
        .par_iter()
        .map(|&(start_voffset, end_voffset)| {
            read_bam_by_interval(bam_path, start_voffset, end_voffset)
        })
        .collect()
}

/// Returns all BAI-derived intervals extended with the tail interval [last_end → EOF].
/// The tail interval captures records beyond the last linear index entry, and is now
/// included in the same parallel pass as the other intervals (Step 5+6).
pub fn get_entire_bam_intervals(
    bam_path: &str,
    intervals: &[(VirtualPosition, VirtualPosition)],
) -> io::Result<Vec<(VirtualPosition, VirtualPosition)>> {
    let &(_, last_end) = intervals
        .last()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "intervals is empty"))?;
    let file_size = File::open(bam_path)?.metadata()?.len();
    let eof_vpos = VirtualPosition::new(file_size, 0)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "BAM file too large"))?;
    let mut all = intervals.to_vec();
    all.push((last_end, eof_vpos));
    Ok(all)
}

pub fn get_entire_bam_reader(
    bam_path: &str,
    intervals: &[(VirtualPosition, VirtualPosition)],
) -> io::Result<Vec<bam::io::Reader<bgzf_io::MultithreadedReader<Cursor<Vec<u8>>>>>> {
    let mut all_interval_readers = read_through_intervals(bam_path, intervals)?;

    // Add the final interval: [last end_coffset, EOF]
    if let Some((_, last_end_voffset)) = intervals.last() {
        let mut bam_file = File::open(bam_path)?;
        bam_file.seek(SeekFrom::Start(last_end_voffset.compressed()))?;
        let mut end_buffer: Vec<u8> = Vec::new();
        bam_file.read_to_end(&mut end_buffer).unwrap_or_else(|e| {
            panic!(
                "read_to_end failed: final voffset {:?}, coffset={}, error={}",
                last_end_voffset,
                last_end_voffset.compressed(),
                e
            )
        });
        let decoder = bgzf_io::MultithreadedReader::with_worker_count(
            NonZero::<usize>::MIN,
            Cursor::new(end_buffer),
        );
        all_interval_readers.push(bam::io::Reader::from(decoder));
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
        timeit(|| read_bam_by_interval(test_bam, intervals[0].0, intervals[0].1))?;
        // println!("{:?}", test_bam);
        Ok(())
    }
}
