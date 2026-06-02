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

/// Counts BAM records whose virtual start position falls in [start, end).
///
/// Unlike `read_bam_by_interval`, this uses a live seeked BGZF reader rather than a
/// pre-loaded byte slice. Pre-loading [start.compressed(), end.compressed()) bytes into a
/// Cursor causes UnexpectedEof whenever a record *starts* before end but its bytes spill
/// into the next BGZF block. With a live reader the block boundary is crossed transparently,
/// and we stop counting by checking virtual_position() *before* each record read.
pub fn count_records_in_virtual_range(
    bam_path: &str,
    start: VirtualPosition,
    end: VirtualPosition,
) -> io::Result<u64> {
    let mut bgzf = bgzf_io::Reader::new(File::open(bam_path)?);
    bgzf.seek(start)?;

    let mut n = 0u64;
    let mut block_size_buf = [0u8; 4];
    let mut record_data = Vec::new();

    loop {
        if bgzf.virtual_position() >= end {
            break;
        }
        // 4-byte LE u32 block_size field that begins every BAM record
        match bgzf.read_exact(&mut block_size_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        }
        let block_size = u32::from_le_bytes(block_size_buf) as usize;
        record_data.resize(block_size, 0);
        bgzf.read_exact(&mut record_data)?;
        n += 1;
    }

    Ok(n)
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
    // Every valid BGZF/BAM file ends with a 28-byte empty EOF block that contains no
    // BAM data. Using file_size as the upper bound includes this block, causing
    // MultithreadedReader to return UnexpectedEof when records() tries to read past it.
    const BGZF_EOF_LEN: u64 = 28;
    let eof_compressed = file_size.saturating_sub(BGZF_EOF_LEN);
    let eof_vpos = VirtualPosition::new(eof_compressed, 0)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "BAM file too large"))?;
    let mut all = intervals.to_vec();
    if last_end.compressed() < eof_compressed {
        all.push((last_end, eof_vpos));
    }
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

    const TEST_BAM: &str = "tests/mt.sorted.bam";
    const TEST_BAI: &str = "tests/mt.sorted.bam.bai";

    // Tail interval is appended and its end sits exactly at file_size - 28 (the BGZF EOF block).
    #[test]
    fn test_get_entire_bam_intervals_appends_tail() -> Result<(), Box<dyn std::error::Error>> {
        let file_size = std::fs::metadata(TEST_BAM)?.len();
        let expected_eof_compressed = file_size - 28;

        let intervals = get_linear_intervals(&get_linear_indexes(TEST_BAI)?)?;
        let all = get_entire_bam_intervals(TEST_BAM, &intervals)?;

        assert_eq!(all.len(), intervals.len() + 1, "tail interval should be appended");
        assert_eq!(
            all.last().unwrap().1.compressed(),
            expected_eof_compressed,
            "tail end must stop before the BGZF EOF block"
        );
        Ok(())
    }

    // When last_end already sits at the BGZF EOF boundary no tail interval is added.
    #[test]
    fn test_get_entire_bam_intervals_skips_tail_at_eof_boundary(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let file_size = std::fs::metadata(TEST_BAM)?.len();
        let eof_compressed = file_size - 28;

        let start = VirtualPosition::new(0, 0).unwrap();
        let last_end = VirtualPosition::new(eof_compressed, 0).unwrap();
        let intervals = vec![(start, last_end)];

        let all = get_entire_bam_intervals(TEST_BAM, &intervals)?;

        assert_eq!(all.len(), intervals.len(), "no tail should be added when last_end is at EOF boundary");
        Ok(())
    }

    // Regression: iterating all intervals must not produce UnexpectedEof.
    // Previously the BGZF EOF block (28 bytes, zero payload) was included in the
    // tail interval, causing MultithreadedReader::records() to fail.
    #[test]
    fn test_get_entire_bam_intervals_no_eof_error() -> Result<(), Box<dyn std::error::Error>> {
        let intervals = get_linear_intervals(&get_linear_indexes(TEST_BAI)?)?;
        let all = get_entire_bam_intervals(TEST_BAM, &intervals)?;

        for (start, end) in all {
            let mut reader = read_bam_by_interval(TEST_BAM, start, end)?;
            for result in reader.records() {
                result?;
            }
        }
        Ok(())
    }

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
