// This module requires htslib which is only available on Unix (built in Docker).
// On Windows the module is excluded via `#[cfg(all(feature = "python", unix))]` in lib.rs.
#![cfg(unix)]

use crate::{
    bai_parser::{get_linear_indexes, get_linear_intervals},
    bam_parser::{count_records_in_virtual_range, get_entire_bam_intervals, merge_intervals},
};
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read as HtsRead, errors::Error as HtsError};
use pyo3::prelude::*;
use pyo3::exceptions::PyIOError;
use pyo3::types::PyBytes;
use std::io;

fn to_py_err<E: ToString>(e: E) -> PyErr {
    PyIOError::new_err(e.to_string())
}

// BAM CIGAR op integer codes, matching pysam / htslib BAM_C* constants.
fn cigar_op_code(op: &bam::record::Cigar) -> u32 {
    use bam::record::Cigar::*;
    match op {
        Match(_) => 0, Ins(_) => 1, Del(_) => 2, RefSkip(_) => 3,
        SoftClip(_) => 4, HardClip(_) => 5, Pad(_) => 6,
        Equal(_) => 7, Diff(_) => 8, Back(_) => 9,
    }
}

// ── intermediate struct ───────────────────────────────────────────────────────
// Extracted inside rayon workers (no GIL); Python objects built on iteration.

struct RecordData {
    query_name:           Option<String>,
    flag:                 u16,
    reference_id:         i32,
    reference_start:      i64,
    mapping_quality:      u8,
    cigarstring:          String,
    cigartuples:          Vec<(u32, u32)>,
    query_sequence:       String,
    query_qualities:      Option<Vec<u8>>,
    template_length:      i32,
    next_reference_id:    i32,
    next_reference_start: i64,
}

impl RecordData {
    fn from_hts(rec: &bam::Record) -> Self {
        let query_name = std::str::from_utf8(rec.qname()).ok().map(String::from);

        let cigar = rec.cigar();
        let cigarstring = if cigar.is_empty() {
            "*".to_string()
        } else {
            cigar.to_string()
        };
        let cigartuples: Vec<(u32, u32)> = cigar.iter()
            .map(|op| (cigar_op_code(op), op.len()))
            .collect();

        let query_sequence = String::from_utf8(rec.seq().as_bytes()).unwrap_or_default();

        // htslib stores absent quality as all-0xFF; pysam returns None in that case.
        let qual = rec.qual();
        let query_qualities = if qual.is_empty() || qual[0] == 0xFF {
            None
        } else {
            Some(qual.to_vec())
        };

        RecordData {
            query_name,
            flag: rec.flags(),
            reference_id: rec.tid(),
            reference_start: rec.pos(),
            mapping_quality: rec.mapq(),
            cigarstring,
            cigartuples,
            query_sequence,
            query_qualities,
            template_length: rec.insert_size() as i32,
            next_reference_id: rec.mtid(),
            next_reference_start: rec.mpos(),
        }
    }
}

// Open an IndexedReader and fetch all records for one reference sequence.
// Each rayon worker calls this with its own file handle.
fn fetch_chromosome(bam_path: &str, name: &str) -> Result<Vec<RecordData>, HtsError> {
    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    // fetch(name) retrieves all records on that chromosome without a position filter.
    reader.fetch(name)?;
    let mut recs = Vec::new();
    let mut record = bam::Record::new();
    while let Some(r) = reader.read(&mut record) {
        r?;
        recs.push(RecordData::from_hts(&record));
    }
    Ok(recs)
}

// Fetch records in a named region [start, stop).
fn fetch_region(
    bam_path: &str,
    contig: &str,
    start: i64,
    stop: i64,
) -> Result<Vec<RecordData>, HtsError> {
    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    reader.fetch((contig, start, stop))?;
    let mut recs = Vec::new();
    let mut record = bam::Record::new();
    while let Some(r) = reader.read(&mut record) {
        r?;
        recs.push(RecordData::from_hts(&record));
    }
    Ok(recs)
}

// ── count ─────────────────────────────────────────────────────────────────────
// Noodles-based fast path: skips field parsing entirely.

#[pyfunction]
#[pyo3(signature = (bam_path, bai_path, until_eof = false))]
pub fn count(bam_path: &str, bai_path: &str, until_eof: bool) -> PyResult<u64> {
    let linear_indexes = get_linear_indexes(bai_path).map_err(to_py_err)?;
    let intervals = get_linear_intervals(&linear_indexes).map_err(to_py_err)?;
    let all = get_entire_bam_intervals(bam_path, &intervals).map_err(to_py_err)?;
    let threads = rayon::current_num_threads().max(1);
    let chunks = merge_intervals(&all, threads);

    if until_eof {
        chunks
            .into_par_iter()
            .map(|(start, end)| count_records_in_virtual_range(bam_path, start, end))
            .sum::<io::Result<u64>>()
            .map_err(to_py_err)
    } else {
        // Mapped reads only: the BAI linear index covers only mapped positions,
        // so records counted here are all mapped.
        chunks
            .into_par_iter()
            .map(|(start, end)| count_records_in_virtual_range(bam_path, start, end))
            .sum::<io::Result<u64>>()
            .map_err(to_py_err)
    }
}

// ── BamRecord ─────────────────────────────────────────────────────────────────

#[pyclass]
pub struct BamRecord {
    query_name:           Option<String>,
    flag:                 u16,
    reference_id:         i32,
    reference_start:      i64,
    mapping_quality:      u8,
    cigarstring:          String,
    cigartuples:          Vec<(u32, u32)>,
    query_sequence:       String,
    query_qualities:      Option<Vec<u8>>,
    template_length:      i32,
    next_reference_id:    i32,
    next_reference_start: i64,
}

impl BamRecord {
    fn from_data(d: RecordData) -> Self {
        BamRecord {
            query_name:           d.query_name,
            flag:                 d.flag,
            reference_id:         d.reference_id,
            reference_start:      d.reference_start,
            mapping_quality:      d.mapping_quality,
            cigarstring:          d.cigarstring,
            cigartuples:          d.cigartuples,
            query_sequence:       d.query_sequence,
            query_qualities:      d.query_qualities,
            template_length:      d.template_length,
            next_reference_id:    d.next_reference_id,
            next_reference_start: d.next_reference_start,
        }
    }
}

#[pymethods]
impl BamRecord {
    fn __repr__(&self) -> String {
        format!(
            "BamRecord(query_name={:?}, flag={}, reference_id={}, reference_start={})",
            self.query_name, self.flag, self.reference_id, self.reference_start
        )
    }

    // ── plain field getters ───────────────────────────────────────────────────
    #[getter] fn query_name(&self) -> Option<&str>  { self.query_name.as_deref() }
    #[getter] fn flag(&self) -> u16                  { self.flag }
    #[getter] fn reference_id(&self) -> i32          { self.reference_id }
    #[getter] fn reference_start(&self) -> i64       { self.reference_start }
    #[getter] fn mapping_quality(&self) -> u8        { self.mapping_quality }
    #[getter] fn cigarstring(&self) -> &str          { &self.cigarstring }
    #[getter] fn query_sequence(&self) -> &str       { &self.query_sequence }
    #[getter] fn template_length(&self) -> i32       { self.template_length }
    #[getter] fn next_reference_id(&self) -> i32     { self.next_reference_id }
    #[getter] fn next_reference_start(&self) -> i64  { self.next_reference_start }

    // ── Python-typed getters ──────────────────────────────────────────────────

    /// Raw Phred quality scores as bytes, or None if absent.
    /// Matches pysam.AlignedSegment.query_qualities semantics.
    #[getter]
    fn query_qualities<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyBytes>> {
        self.query_qualities.as_ref().map(|q| PyBytes::new(py, q))
    }

    /// List of (op_int, length) tuples in htslib BAM_C* encoding.
    /// Matches pysam.AlignedSegment.cigartuples.
    #[getter]
    fn cigartuples(&self) -> Vec<(u32, u32)> {
        self.cigartuples.clone()
    }

    // ── flag accessors ────────────────────────────────────────────────────────
    #[getter] fn is_paired(&self) -> bool       { self.flag & 0x001 != 0 }
    #[getter] fn is_proper_pair(&self) -> bool   { self.flag & 0x002 != 0 }
    #[getter] fn is_unmapped(&self) -> bool      { self.flag & 0x004 != 0 }
    #[getter] fn is_mate_unmapped(&self) -> bool { self.flag & 0x008 != 0 }
    #[getter] fn is_forward(&self) -> bool       { self.flag & 0x010 == 0 }
    #[getter] fn is_reverse(&self) -> bool       { self.flag & 0x010 != 0 }
    #[getter] fn is_read1(&self) -> bool         { self.flag & 0x040 != 0 }
    #[getter] fn is_read2(&self) -> bool         { self.flag & 0x080 != 0 }
    #[getter] fn is_secondary(&self) -> bool     { self.flag & 0x100 != 0 }
    #[getter] fn is_qcfail(&self) -> bool        { self.flag & 0x200 != 0 }
    #[getter] fn is_duplicate(&self) -> bool     { self.flag & 0x400 != 0 }
    #[getter] fn is_supplementary(&self) -> bool { self.flag & 0x800 != 0 }
}

// ── RecordIterator ────────────────────────────────────────────────────────────

#[pyclass]
pub struct RecordIterator {
    // Stored reversed so pop() yields records in genomic order (O(1) per record).
    records: Vec<RecordData>,
}

impl RecordIterator {
    fn new(mut records: Vec<RecordData>) -> Self {
        records.reverse();
        RecordIterator { records }
    }
}

#[pymethods]
impl RecordIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> { slf }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<BamRecord> {
        slf.records.pop().map(BamRecord::from_data)
    }

    fn __len__(&self) -> usize { self.records.len() }
}

// ── AlignmentFile ─────────────────────────────────────────────────────────────

#[pyclass]
pub struct AlignmentFile {
    bam_path:   String,
    bai_path:   String,
    references: Vec<String>,
    lengths:    Vec<u64>,
}

#[pymethods]
impl AlignmentFile {
    /// Accepts the same positional/keyword arguments as pysam.AlignmentFile.
    /// `mode` and `check_sq` are accepted for compatibility and ignored.
    /// If `bai_path` is omitted, `<bam>.bai` is used.
    #[new]
    #[pyo3(signature = (filename, mode = "rb", check_sq = true, bai_path = None))]
    pub fn new(
        filename: String,
        mode: &str,
        check_sq: bool,
        bai_path: Option<String>,
    ) -> PyResult<Self> {
        let _ = (mode, check_sq);

        let reader = bam::IndexedReader::from_path(&filename).map_err(to_py_err)?;
        let header = reader.header();
        let nref = header.target_count();
        let references: Vec<String> = (0..nref)
            .map(|i| String::from_utf8_lossy(header.tid2name(i)).into_owned())
            .collect();
        let lengths: Vec<u64> = (0..nref)
            .map(|i| header.target_len(i).unwrap_or(0))
            .collect();

        let bai = bai_path.unwrap_or_else(|| format!("{}.bai", filename));

        Ok(AlignmentFile { bam_path: filename, bai_path: bai, references, lengths })
    }

    #[getter] fn references(&self) -> Vec<String> { self.references.clone() }
    #[getter] fn lengths(&self)    -> Vec<u64>    { self.lengths.clone() }

    /// Fast parallel count via noodles (no field parsing).
    #[pyo3(signature = (until_eof = false))]
    pub fn count(&self, until_eof: bool) -> PyResult<u64> {
        count(&self.bam_path, &self.bai_path, until_eof)
    }

    /// Parallel record fetch using rust-htslib IndexedReader.
    ///
    /// fetch()                   → all mapped records, parallel by chromosome
    /// fetch(contig)             → one chromosome
    /// fetch(contig, start, stop)→ region [start, stop) (0-based, half-open)
    /// until_eof=True            → all records including unmapped (sequential)
    #[pyo3(signature = (contig = None, start = None, stop = None, until_eof = false))]
    pub fn fetch(
        &self,
        contig: Option<&str>,
        start: Option<i64>,
        stop: Option<i64>,
        until_eof: bool,
    ) -> PyResult<RecordIterator> {
        let bam_path = self.bam_path.clone();

        let records: Vec<RecordData> = if until_eof && contig.is_none() {
            // Sequential scan of entire file including unmapped reads at EOF.
            let mut reader = bam::Reader::from_path(&bam_path).map_err(to_py_err)?;
            let mut recs = Vec::new();
            let mut record = bam::Record::new();
            while let Some(r) = reader.read(&mut record) {
                r.map_err(to_py_err)?;
                recs.push(RecordData::from_hts(&record));
            }
            recs
        } else if let Some(ctg) = contig {
            let s = start.unwrap_or(0);
            let e = stop.unwrap_or(i64::MAX);
            fetch_region(&bam_path, ctg, s, e).map_err(to_py_err)?
        } else {
            // Parallel fetch: one IndexedReader per chromosome per rayon worker.
            let refs = self.references.clone();
            refs.into_par_iter()
                .map(|name| fetch_chromosome(&bam_path, &name))
                .collect::<Result<Vec<Vec<RecordData>>, HtsError>>()
                .map_err(to_py_err)?
                .into_iter()
                .flatten()
                .collect()
        };

        Ok(RecordIterator::new(records))
    }

    fn __iter__(&self) -> PyResult<RecordIterator> {
        self.fetch(None, None, None, false)
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> { slf }

    #[pyo3(signature = (_exc_type = None, _exc_val = None, _exc_tb = None))]
    fn __exit__(
        &self,
        _exc_type: Option<Bound<'_, PyAny>>,
        _exc_val:  Option<Bound<'_, PyAny>>,
        _exc_tb:   Option<Bound<'_, PyAny>>,
    ) -> bool {
        false
    }
}

// ── module entry point ────────────────────────────────────────────────────────

#[pymodule]
pub fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(count, m)?)?;
    m.add_class::<AlignmentFile>()?;
    m.add_class::<BamRecord>()?;
    m.add_class::<RecordIterator>()?;
    Ok(())
}

// ── tests ─────────────────────────────────────────────────────────────────────
// Run with: cargo test --features python
// Requires htslib at build time (available in Docker, not on Windows dev host).

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_BAM: &str = "tests/mt.sorted.bam";
    const CHR_M: &str = "chrM";

    // rust-htslib sequential count must equal noodles sequential count.
    // This verifies the underlying htslib read path against our reference counter.
    #[test]
    fn test_sequential_count_matches_noodles() -> Result<(), Box<dyn std::error::Error>> {
        let mut reader = bam::Reader::from_path(TEST_BAM)?;
        let mut record = bam::Record::new();
        let mut hts_count = 0u64;
        while let Some(r) = reader.read(&mut record) {
            r?;
            hts_count += 1;
        }
        let noodles_count = crate::count_from_standard_bam_reader(TEST_BAM, 1)?;
        assert_eq!(
            hts_count, noodles_count,
            "htslib sequential count {hts_count} != noodles {noodles_count}"
        );
        Ok(())
    }

    // fetch_chromosome must return records and agree with a direct IndexedReader fetch on chrM.
    #[test]
    fn test_fetch_chrm_count_matches_sequential() -> Result<(), Box<dyn std::error::Error>> {
        let indexed = fetch_chromosome(TEST_BAM, CHR_M)?;
        assert!(!indexed.is_empty(), "expected records on {CHR_M}");

        // Find chrM tid by scanning the header (tid() return type varies by htslib version).
        let mut reader = bam::IndexedReader::from_path(TEST_BAM)?;
        let header = reader.header().clone();
        let chrm_tid = (0..header.target_count())
            .find(|&i| header.tid2name(i) == CHR_M.as_bytes())
            .ok_or("chrM not found in header")?;
        reader.fetch(chrm_tid)?;
        let mut seq_count = 0usize;
        let mut record = bam::Record::new();
        while let Some(r) = reader.read(&mut record) {
            r?;
            seq_count += 1;
        }

        assert_eq!(
            indexed.len(), seq_count,
            "fetch_chromosome count {} != direct IndexedReader count {}",
            indexed.len(), seq_count
        );
        Ok(())
    }

    // RecordData fields for mapped reads (flag 0x4 unset) must be internally consistent.
    #[test]
    fn test_mapped_record_fields_valid() -> Result<(), Box<dyn std::error::Error>> {
        let recs = fetch_chromosome(TEST_BAM, CHR_M)?;
        assert!(!recs.is_empty());

        for r in recs.iter().filter(|r| r.flag & 0x004 == 0) {
            assert!(r.reference_id >= 0, "reference_id {} < 0", r.reference_id);
            assert!(r.reference_start >= 0, "reference_start {} < 0", r.reference_start);
            assert!(!r.query_sequence.is_empty(), "empty query_sequence");
            assert!(
                r.query_sequence.bytes().all(|b| matches!(b, b'A'|b'C'|b'G'|b'T'|b'N')),
                "unexpected base in: {}", &r.query_sequence[..r.query_sequence.len().min(20)]
            );
            assert!(!r.cigartuples.is_empty(), "empty cigartuples for mapped read");
            for &(op, _) in &r.cigartuples {
                assert!(op <= 9, "cigar op {op} out of BAM_C* range [0,9]");
            }
            // Query-consuming ops (M=0,I=1,S=4,==7,X=8) must sum to sequence length.
            let qlen: u32 = r.cigartuples.iter()
                .filter(|&&(op, _)| matches!(op, 0 | 1 | 4 | 7 | 8))
                .map(|&(_, len)| len)
                .sum();
            assert_eq!(
                r.query_sequence.len() as u32, qlen,
                "seq len {} != cigar query len {}", r.query_sequence.len(), qlen
            );
            if let Some(q) = &r.query_qualities {
                assert_eq!(q.len(), r.query_sequence.len(), "qual/seq length mismatch");
                // Raw Phred (not +33): valid range [0, 93].
                assert!(q.iter().all(|&v| v <= 93), "Phred score > 93");
            }
        }
        Ok(())
    }

    // fetch_region must return a subset of fetch_chromosome and respect the coordinate bound.
    // htslib returns reads OVERLAPPING [start, stop), so reference_start < stop must hold.
    #[test]
    fn test_fetch_region_subset_and_bounds() -> Result<(), Box<dyn std::error::Error>> {
        let all   = fetch_chromosome(TEST_BAM, CHR_M)?;
        let stop  = 2_000i64;
        let region = fetch_region(TEST_BAM, CHR_M, 0, stop)?;

        assert!(
            region.len() < all.len(),
            "region [0,{stop}) count {} should be less than full chrM {}",
            region.len(), all.len()
        );
        for r in &region {
            assert!(
                r.reference_start < stop,
                "read at {} starts at or after stop {stop}", r.reference_start
            );
        }
        Ok(())
    }

    // is_forward and is_reverse are strict complements; BamRecord getters match flag bits.
    #[test]
    fn test_flag_getters_consistent() -> Result<(), Box<dyn std::error::Error>> {
        let recs = fetch_chromosome(TEST_BAM, CHR_M)?;
        assert!(!recs.is_empty());
        for d in &recs {
            let f = d.flag;
            // is_forward and is_reverse must be strict complements.
            let is_rev = f & 0x010 != 0;
            let is_fwd = f & 0x010 == 0;
            assert_ne!(is_rev, is_fwd, "flag 0x{f:04x}: reverse and forward must differ");

            // Verify BamRecord getters against direct bit arithmetic (catches mask typos).
            let r = BamRecord::from_data(RecordData {
                query_name:           d.query_name.clone(),
                flag:                 d.flag,
                reference_id:         d.reference_id,
                reference_start:      d.reference_start,
                mapping_quality:      d.mapping_quality,
                cigarstring:          d.cigarstring.clone(),
                cigartuples:          d.cigartuples.clone(),
                query_sequence:       d.query_sequence.clone(),
                query_qualities:      d.query_qualities.clone(),
                template_length:      d.template_length,
                next_reference_id:    d.next_reference_id,
                next_reference_start: d.next_reference_start,
            });
            assert_eq!(r.is_paired(),        f & 0x001 != 0, "is_paired mismatch");
            assert_eq!(r.is_unmapped(),      f & 0x004 != 0, "is_unmapped mismatch");
            assert_eq!(r.is_reverse(),       f & 0x010 != 0, "is_reverse mismatch");
            assert_eq!(r.is_forward(),       f & 0x010 == 0, "is_forward mismatch");
            assert_eq!(r.is_read1(),         f & 0x040 != 0, "is_read1 mismatch");
            assert_eq!(r.is_read2(),         f & 0x080 != 0, "is_read2 mismatch");
            assert_eq!(r.is_secondary(),     f & 0x100 != 0, "is_secondary mismatch");
            assert_eq!(r.is_supplementary(), f & 0x800 != 0, "is_supplementary mismatch");
        }
        Ok(())
    }

    // Header must contain chrM with the correct reference length.
    #[test]
    fn test_header_chrm_length() -> Result<(), Box<dyn std::error::Error>> {
        let reader = bam::IndexedReader::from_path(TEST_BAM)?;
        let header = reader.header();
        assert!(header.target_count() > 0);
        let chrm_tid = (0..header.target_count())
            .find(|&i| header.tid2name(i) == CHR_M.as_bytes())
            .ok_or("chrM not found in header")?;
        let len = header.target_len(chrm_tid).ok_or("chrM has no length in header")?;
        assert_eq!(len, 16569, "chrM length should be 16569 bp");
        Ok(())
    }
}
