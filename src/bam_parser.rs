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

/// Merges consecutive intervals into at most `max_chunks` larger intervals.
///
/// Consecutive BAI-derived intervals are always contiguous in compressed space
/// (end_i == start_{i+1}), so merging extends ranges without gaps. Reducing
/// 3000+ intervals to O(thread_count) intervals cuts seek overhead from
/// O(intervals) to O(threads), which is the dominant cost at low thread counts.
///
/// Chunks are sized by compressed byte span so each thread gets roughly equal
/// IO work. Returns the input unchanged when `intervals.len() <= max_chunks`.
pub fn merge_intervals(
    intervals: &[(VirtualPosition, VirtualPosition)],
    max_chunks: usize,
) -> Vec<(VirtualPosition, VirtualPosition)> {
    if intervals.is_empty() || max_chunks == 0 {
        return Vec::new();
    }
    if intervals.len() <= max_chunks {
        return intervals.to_vec();
    }

    let first_start = intervals[0].0.compressed();
    let last_end = intervals[intervals.len() - 1].1.compressed();
    let total_bytes = last_end.saturating_sub(first_start);
    let bytes_per_chunk = (total_bytes + max_chunks as u64 - 1) / max_chunks as u64;

    let mut merged: Vec<(VirtualPosition, VirtualPosition)> = Vec::with_capacity(max_chunks);
    let mut chunk_start = intervals[0].0;
    let mut boundary = first_start + bytes_per_chunk;

    for (i, &(_, end)) in intervals.iter().enumerate() {
        let is_last = i == intervals.len() - 1;
        if end.compressed() >= boundary || is_last {
            merged.push((chunk_start, end));
            if !is_last {
                chunk_start = intervals[i + 1].0;
                boundary = end.compressed() + bytes_per_chunk;
            }
        }
    }

    merged
}

/// Returns the virtual position immediately after the BAM header.
///
/// VirtualPosition(0,0) is a BAI null sentinel that falls inside the BAM header, not
/// on a real record. Any interval whose start is VirtualPosition(0,0) must be adjusted
/// to this position so that count_records_in_virtual_range does not parse header bytes
/// as record data.
pub(crate) fn header_end_vpos(bam_path: &str) -> io::Result<VirtualPosition> {
    let file = File::open(bam_path)?;
    let bgzf = bgzf_io::Reader::new(file);
    let mut bam_reader = bam::io::Reader::from(bgzf);
    bam_reader.read_header()?;
    Ok(bam_reader.get_ref().virtual_position())
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
    // VirtualPosition(0,0) is the BAI null sentinel and always falls within the BAM
    // header. Replace it with the true start of the first record so that
    // count_records_in_virtual_range does not misparse header bytes as a record.
    let zero = VirtualPosition::new(0, 0).unwrap();
    if all.first().map(|&(start, _)| start == zero).unwrap_or(false) {
        all[0].0 = header_end_vpos(bam_path)?;
    }
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

/// Splits a BAM file into parallel read intervals without a BAI index.
///
/// Uses the QuickBAM heuristic (SuppMethods.docx): divide the compressed file
/// into N equal byte ranges, find the BGZF block boundary at each range start
/// via magic-byte scan, decompress that block, and locate the first valid BAM
/// record via `is_valid_bam_record_start`.  Returns VirtualPosition pairs
/// compatible with `count_records_in_virtual_range` and `read_bam_by_interval`.
///
/// `header_end` (from `header_end_vpos`) is the minimum allowed chunk start so
/// that BAM header bytes are never mistaken for record data.
pub fn find_parallel_chunks_no_index(
    bam_path: &str,
    n_chunks: usize,
    header_end: VirtualPosition,
    n_ref: usize,
    ref_lens: &[u32],
) -> io::Result<Vec<(VirtualPosition, VirtualPosition)>> {
    let file_size = File::open(bam_path)?.metadata()?.len() as usize;
    const BGZF_EOF_LEN: usize = 28;
    let data_end = file_size.saturating_sub(BGZF_EOF_LEN);
    let eof_vp = VirtualPosition::new(data_end as u64, 0)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "BAM file too large"))?;

    if n_chunks <= 1 || data_end == 0 {
        return Ok(vec![(header_end, eof_vp)]);
    }

    let bytes_per_chunk = data_end / n_chunks;
    let mut starts: Vec<VirtualPosition> = vec![header_end];

    for i in 1..n_chunks {
        let target = i * bytes_per_chunk;
        if target >= data_end {
            break;
        }

        // Raw compressed bytes: 2× max-block so the found block always fits fully.
        let window_len = (131_072usize + 28).min(data_end - target);
        let mut window = vec![0u8; window_len];
        {
            let mut f = File::open(bam_path)?;
            f.seek(SeekFrom::Start(target as u64))?;
            f.read_exact(&mut window)?;
        }

        let off = match find_next_bgzf_block(&window, 0) {
            Some(p) => p,
            None => continue,
        };
        let block_start = target + off;
        if block_start >= data_end {
            break;
        }

        // BSIZE (total block bytes) and ISIZE (uncompressed bytes) from header/footer.
        if off + 18 > window.len() {
            continue;
        }
        let bsize = u16::from_le_bytes([window[off + 16], window[off + 17]]) as usize + 1;
        if off + bsize > window.len() {
            continue;
        }
        let isize_raw: [u8; 4] = window[off + bsize - 4..off + bsize].try_into().unwrap();
        let isize = u32::from_le_bytes(isize_raw) as usize;
        if isize == 0 {
            continue;
        }

        // Decompress the block via noodles BGZF reader.
        let vp_block = VirtualPosition::new(block_start as u64, 0)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "block offset overflow"))?;
        let mut bgzf = bgzf_io::Reader::new(File::open(bam_path)?);
        bgzf.seek(vp_block)?;
        let mut buf = vec![0u8; isize];
        bgzf.read_exact(&mut buf)?;

        // Scan for the first valid BAM record start.
        let Some(u) = (0..isize).find(|&u| is_valid_bam_record_start(&buf, u, n_ref, ref_lens))
        else {
            continue;
        };
        let Some(vp) = VirtualPosition::new(block_start as u64, u as u16) else {
            continue;
        };

        if vp > header_end && starts.last().copied() != Some(vp) {
            starts.push(vp);
        }
    }

    starts.dedup();

    // Build (start, end) pairs: each chunk ends where the next begins.
    let intervals: Vec<(VirtualPosition, VirtualPosition)> = starts
        .windows(2)
        .map(|w| (w[0], w[1]))
        .chain(starts.last().copied().map(|last| (last, eof_vp)))
        .collect();

    Ok(intervals)
}

/// CIGAR-aware coverage accumulation into a diff array.
///
/// `diff` has length `region_len + 1`; index 0 maps to `region_start` (0-based).
/// M/=/X (ops 0/7/8) consume the reference and are counted as covered.
/// D/N (ops 2/3) skip reference positions — no coverage added.
/// I/S/H/P do not advance ref_pos.
pub(crate) fn apply_cigar_to_diff(
    diff: &mut [i64],
    ref_start: i64,
    cigartuples: &[(u32, u32)],
    region_start: usize,
) {
    let region_end = region_start + diff.len(); // exclusive upper bound
    let mut ref_pos = ref_start;
    for &(op, len) in cigartuples {
        let len = len as i64;
        match op {
            0 | 7 | 8 => {
                let s = ref_pos.max(region_start as i64);
                let e = (ref_pos + len).min(region_end as i64);
                if s < e {
                    let si = (s - region_start as i64) as usize;
                    let ei = (e - region_start as i64) as usize;
                    diff[si] += 1;
                    // ei == diff.len() when the read ends exactly at the region
                    // boundary; the decrement is irrelevant (prefix sum never
                    // reads past the array), so skip rather than panic.
                    if ei < diff.len() {
                        diff[ei] -= 1;
                    }
                }
                ref_pos += len;
            }
            2 | 3 => { ref_pos += len; }
            _ => {}
        }
    }
}

/// Per-position A/C/G/T counts for base_pileup.
///
/// `counts` is a mutable slice of `[u32; 4]` arrays, one per reference position
/// in `[region_start, region_start + counts.len())`. Each array holds
/// `[A_count, C_count, G_count, T_count]`.
///
/// `seq_bytes` is the read sequence as ASCII bytes (from RecordData.query_sequence).
///
/// CIGAR ops M/=/X (0,7,8) consume both ref and query; I/S (1,4) consume query
/// only; D/N (2,3) consume ref only; H/P (5,6) consume neither.
/// Positions outside the region are ignored; reads that extend past the region
/// boundary are correctly clipped.
pub(crate) fn apply_cigar_to_base_counts(
    counts: &mut [[u32; 4]],
    seq_bytes: &[u8],
    ref_start: i64,
    cigartuples: &[(u32, u32)],
    region_start: usize,
) {
    let region_end = region_start + counts.len();
    let mut ref_pos = ref_start;
    let mut query_pos: usize = 0;

    for &(op, len) in cigartuples {
        let len = len as usize;
        match op {
            0 | 7 | 8 => {
                let rp_start = ref_pos as usize;
                let rp_lo = rp_start.max(region_start);
                let rp_hi = (rp_start + len).min(region_end);
                for rp in rp_lo..rp_hi {
                    let qp = query_pos + (rp - rp_start);
                    if let Some(&b) = seq_bytes.get(qp) {
                        let bi: Option<usize> = match b {
                            b'A' | b'a' => Some(0),
                            b'C' | b'c' => Some(1),
                            b'G' | b'g' => Some(2),
                            b'T' | b't' => Some(3),
                            _ => None,
                        };
                        if let Some(bi) = bi {
                            counts[rp - region_start][bi] += 1;
                        }
                    }
                }
                ref_pos += len as i64;
                query_pos += len;
            }
            1 | 4 => { query_pos += len; }
            2 | 3 => { ref_pos += len as i64; }
            5 | 6 => {}
            _ => {}
        }
    }
}

// ── Step 14 building blocks: no-index parallel reading ────────────────────────

/// Scan `data` (raw compressed file bytes) for the next valid BGZF block
/// starting at or after `from`. Checks four non-contiguous magic bytes in the
/// 18-byte BGZF header as described in QuickBAM SuppMethods:
///   byte  0 = 0x1f, byte  1 = 0x8b  (gzip magic)
///   byte 12 = 0x42 ('B'), byte 13 = 0x43 ('C')  (BGZF subfield ID)
/// Returns the compressed offset of the match, or None.
pub fn find_next_bgzf_block(data: &[u8], from: usize) -> Option<usize> {
    let len = data.len();
    let mut i = from;
    while i + 14 <= len {
        if data[i] == 0x1f
            && data[i + 1]  == 0x8b
            && data[i + 12] == 0x42
            && data[i + 13] == 0x43
        {
            return Some(i);
        }
        i += 1;
    }
    None
}

/// Returns whether the bytes at `offset` in a decompressed BAM stream look
/// like a valid BAM record start, using the QuickBAM heuristic criteria
/// (SuppMethods.docx):
///   -1 ≤ ref_id < n_ref
///   ref_id == -1 → pos == -1;  else 0 ≤ pos < ref_lens[ref_id]
///   same constraints for next_ref_id / next_pos
///   read_name[l_read_name - 1] == 0  (null-terminated)
///
/// BAM record layout at `offset` (all LE):
///   [0..4]   block_size  u32   (does NOT include itself)
///   [4..8]   ref_id      i32
///   [8..12]  pos         i32   (0-based; -1 for unmapped)
///   [12]     l_read_name u8    (includes NUL terminator)
///   [13]     mapq        u8
///   [14..16] bin         u16
///   [16..18] n_cigar_op  u16
///   [18..20] flag        u16
///   [20..24] l_seq       i32
///   [24..28] next_ref    i32
///   [28..32] next_pos    i32
///   [32..36] tlen        i32
///   [36..]   read_name (l_read_name bytes, NUL-terminated)
pub fn is_valid_bam_record_start(
    data: &[u8],
    offset: usize,
    n_ref: usize,
    ref_lens: &[u32],
) -> bool {
    if offset + 36 > data.len() {
        return false;
    }
    let d = &data[offset..];

    let ref_id      = i32::from_le_bytes(d[4..8].try_into().unwrap());
    let pos         = i32::from_le_bytes(d[8..12].try_into().unwrap());
    let l_read_name = d[12] as usize;
    let next_ref    = i32::from_le_bytes(d[24..28].try_into().unwrap());
    let next_pos    = i32::from_le_bytes(d[28..32].try_into().unwrap());

    // ref_id bounds
    if ref_id < -1 || (ref_id >= 0 && ref_id as usize >= n_ref) {
        return false;
    }
    // pos: unmapped sentinel is -1; mapped reads must be in range
    if ref_id == -1 {
        if pos != -1 { return false; }
    } else if pos < 0 || pos as u32 >= ref_lens[ref_id as usize] {
        return false;
    }
    // next_ref bounds
    if next_ref < -1 || (next_ref >= 0 && next_ref as usize >= n_ref) {
        return false;
    }
    // next_pos: same rules
    if next_ref == -1 {
        if next_pos != -1 { return false; }
    } else if next_pos < 0 || next_pos as u32 >= ref_lens[next_ref as usize] {
        return false;
    }
    // read_name must be non-empty and null-terminated
    if l_read_name == 0 || offset + 36 + l_read_name > data.len() {
        return false;
    }
    data[offset + 36 + l_read_name - 1] == 0
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

    // count_records_in_virtual_range must agree with the standard sequential reader.
    // This is the primary correctness gate for bench_count: if it fails here the
    // parallel count against full.bam will also be wrong.
    #[test]
    fn test_count_records_in_virtual_range_matches_standard(
    ) -> Result<(), Box<dyn std::error::Error>> {
        use rayon::prelude::*;

        let linear_indexes = get_linear_indexes(TEST_BAI)?;

        // Verify: are there VirtualPosition::MIN (zero) entries in the raw linear index?
        // Zero entries (empty BAI windows) must not be fed to count_records_in_virtual_range
        // as the first interval, or the BAM header bytes get parsed as a fake record.
        let has_zero = linear_indexes
            .iter()
            .any(|vp| *vp == VirtualPosition::new(0, 0).unwrap());
        // Print for diagnostic purposes (does not fail the test).
        println!("linear index contains VirtualPosition(0,0) entries: {has_zero}");

        let intervals = get_linear_intervals(&linear_indexes)?;
        let all = get_entire_bam_intervals(TEST_BAM, &intervals)?;

        let parallel_count: u64 = all
            .into_par_iter()
            .map(|(start, end)| count_records_in_virtual_range(TEST_BAM, start, end))
            .sum::<io::Result<u64>>()?;

        let standard_count = crate::count_from_standard_bam_reader(TEST_BAM, 1)?;
        assert_eq!(
            parallel_count, standard_count,
            "count_records_in_virtual_range total must match standard reader"
        );
        Ok(())
    }

    #[test]
    fn test_merge_intervals_reduces_to_max_chunks() {
        let vp = |c: u64| VirtualPosition::new(c, 0).unwrap();
        // 10 contiguous intervals spanning compressed offsets 0..10
        let intervals: Vec<(VirtualPosition, VirtualPosition)> =
            (0u64..10).map(|i| (vp(i), vp(i + 1))).collect();

        let merged = merge_intervals(&intervals, 3);
        assert!(merged.len() <= 3, "should produce at most 3 chunks");
        // First chunk starts at same position as first interval
        assert_eq!(merged[0].0.compressed(), 0);
        // Last chunk ends at same position as last interval
        assert_eq!(merged.last().unwrap().1.compressed(), 10);
    }

    #[test]
    fn test_merge_intervals_noop_when_fewer_than_max() {
        let vp = |c: u64| VirtualPosition::new(c, 0).unwrap();
        let intervals: Vec<_> = (0u64..5).map(|i| (vp(i), vp(i + 1))).collect();
        let merged = merge_intervals(&intervals, 10);
        assert_eq!(merged.len(), intervals.len());
    }

    #[test]
    fn test_merge_intervals_count_matches_standard() -> Result<(), Box<dyn std::error::Error>> {
        use rayon::prelude::*;

        let linear_indexes = get_linear_indexes(TEST_BAI)?;
        let intervals = get_linear_intervals(&linear_indexes)?;
        let all = get_entire_bam_intervals(TEST_BAM, &intervals)?;

        // Merge into 2 chunks (stress-tests large merged intervals)
        let merged = merge_intervals(&all, 2);
        assert!(merged.len() <= 2);

        let merged_count: u64 = merged
            .into_par_iter()
            .map(|(start, end)| count_records_in_virtual_range(TEST_BAM, start, end))
            .sum::<io::Result<u64>>()?;

        let standard_count = crate::count_from_standard_bam_reader(TEST_BAM, 1)?;
        assert_eq!(
            merged_count, standard_count,
            "merged intervals count must match standard reader"
        );
        Ok(())
    }

    // ── Step 14: no-index parallel count ──────────────────────────────────────

    #[test]
    fn test_count_no_index_matches_sequential() -> Result<(), Box<dyn std::error::Error>> {
        use rayon::prelude::*;

        let header = crate::get_bam_header(TEST_BAM)?;
        let n_ref = header.reference_sequences().len();
        let ref_lens: Vec<u32> = header
            .reference_sequences()
            .values()
            .map(|rs| rs.length().get() as u32)
            .collect();
        let hdr_end = header_end_vpos(TEST_BAM)?;

        let chunks = find_parallel_chunks_no_index(TEST_BAM, 4, hdr_end, n_ref, &ref_lens)?;
        assert!(!chunks.is_empty(), "should produce at least one chunk");

        let count_no_idx: u64 = chunks
            .into_par_iter()
            .map(|(start, end)| count_records_in_virtual_range(TEST_BAM, start, end))
            .sum::<io::Result<u64>>()?;

        let sequential = crate::count_from_standard_bam_reader(TEST_BAM, 1)?;
        assert_eq!(
            count_no_idx, sequential,
            "no-index parallel count {count_no_idx} != sequential {sequential}"
        );
        Ok(())
    }

    // ── coverage: apply_cigar_to_diff ─────────────────────────────────────────

    // Simple read: "100M" starting at pos 0 → coverage[0..100] = 1, rest = 0.
    #[test]
    fn test_apply_cigar_100m_full_coverage() {
        let mut diff = vec![0i64; 201]; // region [0, 200), +1 sentinel
        apply_cigar_to_diff(&mut diff, 0, &[(0, 100)], 0);
        let mut cov = vec![0u32; 200];
        let mut run = 0i64;
        for (i, &d) in diff[..200].iter().enumerate() {
            run += d;
            cov[i] = run.max(0) as u32;
        }
        assert!(cov[..100].iter().all(|&d| d == 1), "positions 0..100 must be covered");
        assert!(cov[100..].iter().all(|&d| d == 0), "positions 100..200 must be 0");
    }

    // Read ending exactly at region boundary must not panic (off-by-one guard).
    #[test]
    fn test_apply_cigar_read_ends_at_region_boundary() {
        // Region [0, 10), diff len = 10. A "10M" read starting at 0 ends at 10
        // which equals region_end — the decrement falls on index 10 (== len).
        let mut diff = vec![0i64; 10];
        apply_cigar_to_diff(&mut diff, 0, &[(0, 10)], 0); // must not panic
        let cov: Vec<u32> = diff.iter().scan(0i64, |run, &d| { *run += d; Some((*run).max(0) as u32) }).collect();
        assert!(cov.iter().all(|&d| d == 1), "all 10 positions must be covered");
    }

    // "50M2D50M" — deletion is NOT counted as covered.
    #[test]
    fn test_apply_cigar_deletion_skipped() {
        let mut diff = vec![0i64; 153]; // region [0, 152), +1 sentinel
        apply_cigar_to_diff(&mut diff, 0, &[(0, 50), (2, 2), (0, 50)], 0);
        let mut cov = vec![0u32; 152];
        let mut run = 0i64;
        for (i, &d) in diff[..152].iter().enumerate() {
            run += d;
            cov[i] = run.max(0) as u32;
        }
        assert!(cov[..50].iter().all(|&d| d == 1),   "positions 0..50 must be covered");
        assert!(cov[50..52].iter().all(|&d| d == 0),  "deletion 50..52 must NOT be covered");
        assert!(cov[52..102].iter().all(|&d| d == 1), "positions 52..102 must be covered");
        assert!(cov[102..].iter().all(|&d| d == 0),   "tail must be 0");
    }

    // base_counts: "10M" read "ACGTACGTAC" → each position has exactly one base.
    #[test]
    fn test_apply_cigar_to_base_counts_simple() {
        let seq = b"ACGTACGTAC";
        let mut counts = vec![[0u32; 4]; 10];
        apply_cigar_to_base_counts(&mut counts, seq, 0, &[(0, 10)], 0);
        let expected = [
            [1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1],
            [1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1],
            [1,0,0,0], [0,1,0,0],
        ];
        assert_eq!(counts, expected);
        assert_eq!(counts.iter().map(|c| c.iter().sum::<u32>()).sum::<u32>(), 10);
    }

    // base_counts: deletion skipped — D does not produce base counts.
    #[test]
    fn test_apply_cigar_to_base_counts_deletion_skipped() {
        // "5M2D5M" — reads 5 bases, skips 2 ref, reads 5 more bases.
        let seq = b"AAAAATTTTT";
        let mut counts = vec![[0u32; 4]; 12];
        apply_cigar_to_base_counts(&mut counts, seq, 0, &[(0, 5), (2, 2), (0, 5)], 0);
        let total: u32 = counts.iter().map(|c| c.iter().sum::<u32>()).sum();
        assert_eq!(total, 10);             // 10 query bases counted
        assert_eq!(counts[5], [0,0,0,0]); // deletion pos 5 has no bases
        assert_eq!(counts[6], [0,0,0,0]); // deletion pos 6 has no bases
    }

    // base_counts: region clipping — read extending past region boundary is ignored.
    #[test]
    fn test_apply_cigar_to_base_counts_boundary() {
        // Read at pos 8, "5M", region [0, 10). Bases at 8,9 count; 10..12 ignored.
        let seq = b"ACGTA";
        let mut counts = vec![[0u32; 4]; 10];
        apply_cigar_to_base_counts(&mut counts, seq, 8, &[(0, 5)], 0);
        let total: u32 = counts.iter().map(|c| c.iter().sum::<u32>()).sum();
        assert_eq!(total, 2);             // only positions 8,9 in region
        assert_eq!(counts[8], [1,0,0,0]); // A
        assert_eq!(counts[9], [0,1,0,0]); // C
    }

    // Reads with real BAM: Σcov == Σ(M/X/= bases) for all mapped reads on chrM.
    #[test]
    fn test_coverage_sum_equals_aligned_bases() -> Result<(), Box<dyn std::error::Error>> {
        use noodles::bam as nb;
        use noodles::core::Region;
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::alignment::record::Cigar as CigarTrait;
        use noodles::sam::alignment::record::cigar::op::Kind;

        let header = crate::get_bam_header(TEST_BAM)?;
        let contig_len = header
            .reference_sequences()
            .get(b"chrM".as_ref())
            .map(|rs| rs.length().get())
            .unwrap();

        let region: Region = "chrM".parse().unwrap();
        let index = nb::bai::fs::read(TEST_BAI)?;
        let mut reader = nb::io::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_path(TEST_BAM)?;
        let _ = reader.read_header()?;

        let mut diff = vec![0i64; contig_len + 1];
        let mut aligned_bases: u64 = 0;

        for result in reader.query(&header, &region)?.records() {
            let rec = result?;
            let buf = RecordBuf::try_from_alignment_record(&header, &rec)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            if buf.flags().is_unmapped() { continue; }
            let ref_start = match buf.alignment_start() {
                Some(pos) => usize::from(pos) as i64 - 1,
                None => continue,
            };
            let tuples: Vec<(u32, u32)> = buf.cigar().iter()
                .filter_map(|r| r.ok())
                .map(|op| {
                    let code = match op.kind() {
                        Kind::Match            => 0,
                        Kind::Insertion        => 1,
                        Kind::Deletion         => 2,
                        Kind::Skip             => 3,
                        Kind::SoftClip         => 4,
                        Kind::HardClip         => 5,
                        Kind::Pad              => 6,
                        Kind::SequenceMatch    => 7,
                        Kind::SequenceMismatch => 8,
                    };
                    (code, op.len() as u32)
                })
                .collect();
            aligned_bases += tuples.iter()
                .filter(|&&(op, _)| matches!(op, 0 | 7 | 8))
                .map(|&(_, l)| l as u64)
                .sum::<u64>();
            apply_cigar_to_diff(&mut diff, ref_start, &tuples, 0);
        }

        let cov_sum: u64 = {
            let mut run = 0i64;
            (0..contig_len).map(|i| { run += diff[i]; run.max(0) as u64 }).sum()
        };
        assert_eq!(
            cov_sum, aligned_bases,
            "Σcov ({cov_sum}) != Σaligned_bases ({aligned_bases})"
        );
        assert!(cov_sum > 0, "expected non-zero coverage on chrM");
        Ok(())
    }

    // ── Step 14 sanity checks ─────────────────────────────────────────────────

    // Verify find_next_bgzf_block can traverse the entire file by jumping each
    // block via its BSIZE field — a full chain walk proves the scanner finds
    // real block boundaries, not false positives.
    #[test]
    fn test_bgzf_scan_covers_whole_file() -> Result<(), Box<dyn std::error::Error>> {
        let compressed = std::fs::read(TEST_BAM)?;
        let file_size  = compressed.len();

        let mut i     = 0usize;
        let mut count = 0usize;
        while i < file_size {
            let pos = find_next_bgzf_block(&compressed, i)
                .ok_or("no BGZF block found")?;
            assert_eq!(pos, i, "scanner jumped past a block at offset {i}");
            // BSIZE field is at bytes [16..18] of the header (LE u16).
            let bsize = u16::from_le_bytes([compressed[pos + 16], compressed[pos + 17]]) as usize + 1;
            i = pos + bsize;
            count += 1;
        }
        assert!(count > 5, "expected >5 BGZF blocks, got {count}");
        Ok(())
    }

    // Verify is_valid_bam_record_start accepts all known-good record starts
    // and has a negligible false-positive rate at non-record offsets.
    #[test]
    fn test_heuristic_accepts_real_records() -> Result<(), Box<dyn std::error::Error>> {
        use std::io::Read;

        // Decompress the whole BAM to raw bytes.
        let file = File::open(TEST_BAM)?;
        let mut bgzf = noodles::bgzf::io::Reader::new(file);
        let mut raw = Vec::new();
        bgzf.read_to_end(&mut raw)?;

        // Parse BAM header to find where records start and get ref info.
        let header = crate::get_bam_header(TEST_BAM)?;
        let n_ref: usize = header.reference_sequences().len();
        let ref_lens: Vec<u32> = header.reference_sequences().values()
            .map(|rs| rs.length().get() as u32)
            .collect();

        // Walk the raw BAM header manually to find byte offset of first record.
        // Layout: magic(4) + l_text(4) + text(l_text) + n_ref(4) + [l_name(4)+name+l_ref(4)]×n
        let mut off = 4usize; // skip "BAM\1"
        let l_text = u32::from_le_bytes(raw[off..off+4].try_into()?) as usize;
        off += 4 + l_text;
        let n_ref_hdr = u32::from_le_bytes(raw[off..off+4].try_into()?) as usize;
        off += 4;
        for _ in 0..n_ref_hdr {
            let l_name = u32::from_le_bytes(raw[off..off+4].try_into()?) as usize;
            off += 4 + l_name + 4;
        }
        // `off` now points to the first BAM record.

        // 1. Every known-good record start must pass the heuristic.
        let mut rec_off = off;
        let mut checked = 0usize;
        while checked < 20 && rec_off + 4 <= raw.len() {
            assert!(
                is_valid_bam_record_start(&raw, rec_off, n_ref, &ref_lens),
                "heuristic rejected known-good record at decompressed offset {rec_off}"
            );
            let block_size = u32::from_le_bytes(raw[rec_off..rec_off+4].try_into()?) as usize;
            rec_off += 4 + block_size;
            checked += 1;
        }
        assert!(checked > 0, "no records checked");

        // 2. False-positive rate at non-record positions should be very low.
        let window = (off + 1)..(off + 200).min(raw.len());
        let false_positives = window
            .filter(|&p| is_valid_bam_record_start(&raw, p, n_ref, &ref_lens))
            .count();
        assert!(
            false_positives < 3,
            "too many false positives in first record body: {false_positives}"
        );

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
