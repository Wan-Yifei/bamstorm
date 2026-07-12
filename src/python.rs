// This module requires htslib which is only available on Unix (built in Docker).
// On Windows the module is excluded via `#[cfg(all(feature = "python", unix))]` in lib.rs.
#![cfg(unix)]

use crate::{
    bai_parser::{get_linear_indexes, get_linear_intervals},
    bam_parser::{get_entire_bam_intervals, merge_intervals, read_bam_by_interval},
};
use noodles::bam as noodles_bam;
use noodles::core::Region;
use noodles::sam as noodles_sam;
use noodles::sam::alignment::RecordBuf;
use rayon::prelude::*;
use pyo3::prelude::*;
use pyo3::exceptions::{PyIOError, PyKeyError, PyValueError};
use pyo3::types::PyBytes;
use std::fs::File;
use std::io;

fn to_py_err<E: ToString>(e: E) -> PyErr {
    PyIOError::new_err(e.to_string())
}

// ── tag value storage ─────────────────────────────────────────────────────────
// Owned representation of a BAM auxiliary tag value.

#[derive(Clone, Debug)]
enum TagOwned {
    Int(i64),
    Float(f32),
    Str(String),
    IntArray(Vec<i64>),
    FloatArray(Vec<f32>),
}

fn tag_owned_to_py(py: Python<'_>, v: &TagOwned) -> PyObject {
    match v {
        TagOwned::Int(i)         => i.into_py(py),
        TagOwned::Float(f)       => f.into_py(py),
        TagOwned::Str(s)         => s.into_py(py),
        TagOwned::IntArray(a)    => a.clone().into_py(py),
        TagOwned::FloatArray(a)  => a.clone().into_py(py),
    }
}

// ── CIGAR helpers ─────────────────────────────────────────────────────────────

fn cigar_kind_to_code(kind: noodles_sam::alignment::record::cigar::op::Kind) -> u32 {
    use noodles_sam::alignment::record::cigar::op::Kind::*;
    match kind {
        Match => 0, Insertion => 1, Deletion => 2, Skip => 3,
        SoftClip => 4, HardClip => 5, Pad => 6,
        SequenceMatch => 7, SequenceMismatch => 8,
    }
}

const CIGAR_CHARS: [char; 9] = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'];

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
    tags:                 Vec<([u8; 2], TagOwned)>,
}

impl RecordData {
    // Noodles path: used for all fetch paths.
    // RecordBuf converts the lazy BAM record into fully-owned fields including tags.
    fn from_noodles(
        rec: &noodles_bam::Record,
        header: &noodles_sam::Header,
    ) -> io::Result<Self> {
        // Bring SAM alignment record traits into scope so their methods resolve.
        use noodles_sam::alignment::record::Cigar as CigarTrait;
        use noodles_sam::alignment::record::Sequence as SequenceTrait;
        // The record_buf's Data uses record_buf::data::field::Value, distinct from
        // record::data::field::Value.
        use noodles_sam::alignment::record_buf::data::field::Value as TagValue;

        let buf = RecordBuf::try_from_alignment_record(header, rec)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let query_name = buf.name()
            .and_then(|n| std::str::from_utf8(n.as_ref()).ok())
            .map(String::from);

        let flag = buf.flags().bits();

        let reference_id = buf.reference_sequence_id()
            .map(|id| id as i32)
            .unwrap_or(-1);

        let reference_start = buf.alignment_start()
            .map(|pos| usize::from(pos) as i64 - 1)
            .unwrap_or(-1);

        let mapping_quality = buf.mapping_quality()
            .map(u8::from)
            .unwrap_or(255);

        // CigarTrait::iter() yields io::Result<Op>; filter_map drops any decode errors.
        let cigartuples: Vec<(u32, u32)> = buf.cigar()
            .iter()
            .filter_map(|r| r.ok())
            .map(|op| (cigar_kind_to_code(op.kind()), op.len() as u32))
            .collect();

        let cigarstring = if cigartuples.is_empty() {
            "*".to_string()
        } else {
            cigartuples.iter()
                .map(|&(op, len)| format!("{}{}", len, CIGAR_CHARS[op as usize]))
                .collect()
        };

        // SequenceTrait has no iter(); use indexed get().
        let seq = buf.sequence();
        let query_sequence: String = (0..seq.len())
            .filter_map(|i| seq.get(i))
            .map(char::from)
            .collect();

        // record_buf::QualityScores::iter() yields u8 by value.
        let qual: Vec<u8> = buf.quality_scores().iter().collect();
        let query_qualities = if qual.is_empty() || qual.first() == Some(&0xFF) {
            None
        } else {
            Some(qual)
        };

        let template_length = buf.template_length();

        let next_reference_id = buf.mate_reference_sequence_id()
            .map(|id| id as i32)
            .unwrap_or(-1);

        let next_reference_start = buf.mate_alignment_start()
            .map(|pos| usize::from(pos) as i64 - 1)
            .unwrap_or(-1);

        let tags = buf.data().iter()
            .map(|(tag, value)| {
                // as_ref() gives &[u8; 2]; dereference to copy the 2-byte array.
                let tag_bytes: [u8; 2] = *tag.as_ref();

                let owned = match value {
                    TagValue::Character(c) => TagOwned::Str(c.to_string()),
                    TagValue::Int8(v)   => TagOwned::Int(*v as i64),
                    TagValue::UInt8(v)  => TagOwned::Int(*v as i64),
                    TagValue::Int16(v)  => TagOwned::Int(*v as i64),
                    TagValue::UInt16(v) => TagOwned::Int(*v as i64),
                    TagValue::Int32(v)  => TagOwned::Int(*v as i64),
                    TagValue::UInt32(v) => TagOwned::Int(*v as i64),
                    TagValue::Float(v)  => TagOwned::Float(*v),
                    TagValue::String(s) => TagOwned::Str(s.to_string()),
                    TagValue::Hex(s)    => TagOwned::Str(s.to_string()),
                    TagValue::Array(arr) => {
                        use noodles_sam::alignment::record_buf::data::field::value::Array;
                        // record_buf Array variants hold owned plain values.
                        match arr {
                            Array::Int8(vs)   => TagOwned::IntArray(vs.iter().map(|&v| v as i64).collect()),
                            Array::UInt8(vs)  => TagOwned::IntArray(vs.iter().map(|&v| v as i64).collect()),
                            Array::Int16(vs)  => TagOwned::IntArray(vs.iter().map(|&v| v as i64).collect()),
                            Array::UInt16(vs) => TagOwned::IntArray(vs.iter().map(|&v| v as i64).collect()),
                            Array::Int32(vs)  => TagOwned::IntArray(vs.iter().map(|&v| v as i64).collect()),
                            Array::UInt32(vs) => TagOwned::IntArray(vs.iter().map(|&v| v as i64).collect()),
                            Array::Float(vs)  => TagOwned::FloatArray(vs.iter().copied().collect()),
                        }
                    }
                };
                (tag_bytes, owned)
            })
            .collect();

        Ok(RecordData {
            query_name,
            flag,
            reference_id,
            reference_start,
            mapping_quality,
            cigarstring,
            cigartuples,
            query_sequence,
            query_qualities,
            template_length,
            next_reference_id,
            next_reference_start,
            tags,
        })
    }
}

// Fetch records overlapping a contig or region via noodles indexed reader.
// Coordinates follow pysam convention: 0-based half-open [start, stop).
// When both start and stop are None the entire contig is returned.
fn fetch_contig_or_region(
    bam_path: &str,
    bai_path: &str,
    header: &noodles_sam::Header,
    contig: &str,
    start: Option<i64>,
    stop: Option<i64>,
) -> io::Result<Vec<RecordData>> {
    // Build noodles Region (1-based closed intervals).
    // pysam [s, e) 0-based → noodles [s+1, e] 1-based.
    let region: Region = match (start, stop) {
        (None, None) => contig.parse()
            .map_err(|e: noodles::core::region::ParseError|
                io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?,
        (s, e) => {
            let s1 = s.unwrap_or(0) + 1;
            let region_str = match e {
                Some(end) => format!("{contig}:{s1}-{end}"),
                None      => format!("{contig}:{s1}"),
            };
            region_str.parse()
                .map_err(|e: noodles::core::region::ParseError|
                    io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?
        }
    };

    let index = noodles_bam::bai::fs::read(bai_path)?;
    let mut reader = noodles_bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(bam_path)?;
    let _ = reader.read_header()?;

    let mut recs = Vec::new();
    for result in reader.query(header, &region)?.records() {
        recs.push(RecordData::from_noodles(&result?, header)?);
    }
    Ok(recs)
}

// ── count ─────────────────────────────────────────────────────────────────────

#[pyfunction]
#[pyo3(signature = (bam_path, bai_path, until_eof = false))]
pub fn count(bam_path: &str, bai_path: &str, until_eof: bool) -> PyResult<u64> {
    use crate::bam_parser::count_records_in_virtual_range;
    let linear_indexes = get_linear_indexes(bai_path).map_err(to_py_err)?;
    let intervals = get_linear_intervals(&linear_indexes).map_err(to_py_err)?;
    let all = get_entire_bam_intervals(bam_path, &intervals).map_err(to_py_err)?;
    let threads = rayon::current_num_threads().max(1);
    let chunks = merge_intervals(&all, threads);

    let _ = until_eof;
    chunks
        .into_par_iter()
        .map(|(start, end)| count_records_in_virtual_range(bam_path, start, end))
        .sum::<io::Result<u64>>()
        .map_err(to_py_err)
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
    tags:                 Vec<([u8; 2], TagOwned)>,
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
            tags:                 d.tags,
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
    #[getter]
    fn query_qualities<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyBytes>> {
        self.query_qualities.as_ref().map(|q| PyBytes::new_bound(py, q))
    }

    /// List of (op_int, length) tuples in htslib BAM_C* encoding.
    #[getter]
    fn cigartuples(&self) -> Vec<(u32, u32)> {
        self.cigartuples.clone()
    }

    // ── tag accessors ─────────────────────────────────────────────────────────

    /// Return the value of an auxiliary tag by its 2-character name.
    /// Raises KeyError if the tag is absent, ValueError if the name is not 2 chars.
    fn get_tag(&self, py: Python<'_>, tag: &str) -> PyResult<PyObject> {
        let tag_bytes: [u8; 2] = tag.as_bytes()
            .try_into()
            .map_err(|_| PyValueError::new_err("tag name must be exactly 2 characters"))?;
        for (t, v) in &self.tags {
            if *t == tag_bytes {
                return Ok(tag_owned_to_py(py, v));
            }
        }
        Err(PyKeyError::new_err(format!("tag '{}' not found", tag)))
    }

    /// All auxiliary tags as a list of (name, value) tuples.
    #[getter]
    fn tags(&self, py: Python<'_>) -> Vec<(String, PyObject)> {
        self.tags.iter()
            .map(|(t, v)| {
                let key = String::from_utf8_lossy(t).into_owned();
                (key, tag_owned_to_py(py, v))
            })
            .collect()
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
    #[new]
    #[pyo3(signature = (filename, mode = "rb", check_sq = true, bai_path = None))]
    pub fn new(
        filename: String,
        mode: &str,
        check_sq: bool,
        bai_path: Option<String>,
    ) -> PyResult<Self> {
        let _ = (mode, check_sq);

        // Use noodles header reader — no rust-htslib dependency here.
        let header = crate::get_bam_header(&filename).map_err(to_py_err)?;
        let references: Vec<String> = header.reference_sequences().keys()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .collect();
        let lengths: Vec<u64> = header.reference_sequences().values()
            .map(|rs| rs.length().get() as u64)
            .collect();

        let bai = bai_path.unwrap_or_else(|| format!("{}.bai", filename));
        Ok(AlignmentFile { bam_path: filename, bai_path: bai, references, lengths })
    }

    #[getter] fn references(&self) -> Vec<String> { self.references.clone() }
    #[getter] fn lengths(&self)    -> Vec<u64>    { self.lengths.clone() }

    #[pyo3(signature = (until_eof = false))]
    pub fn count(&self, until_eof: bool) -> PyResult<u64> {
        count(&self.bam_path, &self.bai_path, until_eof)
    }

    /// Parallel record fetch.
    ///
    /// fetch()                    → all mapped records via BAI intervals (noodles, tags included)
    /// fetch(contig)              → one chromosome via BAI bin-index (noodles, tags included)
    /// fetch(contig, start, stop) → region [start, stop) 0-based (noodles, tags included)
    /// until_eof=True             → all records including unmapped, sequential (noodles, tags included)
    #[pyo3(signature = (contig = None, start = None, stop = None, until_eof = false))]
    pub fn fetch(
        &self,
        contig: Option<&str>,
        start: Option<i64>,
        stop: Option<i64>,
        until_eof: bool,
    ) -> PyResult<RecordIterator> {
        let bam_path = self.bam_path.clone();
        let bai_path = self.bai_path.clone();

        let records: Vec<RecordData> = if let Some(ctg) = contig {
            // Region / chromosome fetch via noodles indexed reader (BAI bin-index query).
            let header = crate::get_bam_header(&bam_path).map_err(to_py_err)?;
            fetch_contig_or_region(&bam_path, &bai_path, &header, ctg, start, stop)
                .map_err(to_py_err)?
        } else if until_eof {
            // Sequential noodles reader — includes unmapped reads at EOF.
            let file = File::open(&bam_path).map_err(to_py_err)?;
            let mut reader = noodles_bam::io::Reader::new(file);
            let header = reader.read_header().map_err(to_py_err)?;
            let mut recs = Vec::new();
            for result in reader.records() {
                let rec = result.map_err(to_py_err)?;
                recs.push(RecordData::from_noodles(&rec, &header).map_err(to_py_err)?);
            }
            recs
        } else {
            // Parallel fetch via BAI linear intervals + noodles BGZF reader.
            // This is the high-performance path (same parallelism as count()).
            let linear_indexes = get_linear_indexes(&bai_path).map_err(to_py_err)?;
            let intervals = get_linear_intervals(&linear_indexes).map_err(to_py_err)?;
            let all = get_entire_bam_intervals(&bam_path, &intervals).map_err(to_py_err)?;
            let threads = rayon::current_num_threads().max(1);
            let chunks = merge_intervals(&all, threads);
            let header = crate::get_bam_header(&bam_path).map_err(to_py_err)?;

            chunks
                .into_par_iter()
                .map(|(start_vp, end_vp)| -> io::Result<Vec<RecordData>> {
                    let mut reader = read_bam_by_interval(&bam_path, start_vp, end_vp)?;
                    let mut recs = Vec::new();
                    for result in reader.records() {
                        recs.push(RecordData::from_noodles(&result?, &header)?);
                    }
                    Ok(recs)
                })
                .collect::<io::Result<Vec<Vec<RecordData>>>>()
                .map_err(to_py_err)?
                .into_iter()
                .flatten()
                .collect()
        };

        Ok(RecordIterator::new(records))
    }

    /// Returns per-reference mapping statistics read directly from the BAI index.
    /// Output: list of (name, length, n_mapped, n_unmapped).
    /// Last entry is ("*", 0, 0, unplaced_unmapped) — reads with no reference.
    /// Equivalent to `samtools idxstats`, O(1) — does not read the BAM file.
    pub fn idxstats(&self) -> PyResult<Vec<(String, u64, u64, u64)>> {
        use noodles::csi::BinningIndex;
        use noodles::csi::binning_index::ReferenceSequence as _;

        let index = noodles_bam::bai::fs::read(&self.bai_path).map_err(to_py_err)?;

        let mut result: Vec<(String, u64, u64, u64)> = self.references.iter()
            .zip(self.lengths.iter())
            .zip(index.reference_sequences())
            .map(|((name, &len), ref_seq)| {
                let (n_mapped, n_unmapped) = ref_seq.metadata()
                    .map(|m| (m.mapped_record_count(), m.unmapped_record_count()))
                    .unwrap_or((0, 0));
                (name.clone(), len, n_mapped, n_unmapped)
            })
            .collect();

        let unplaced = index.unplaced_unmapped_record_count().unwrap_or(0);
        result.push(("*".to_string(), 0, 0, unplaced));

        Ok(result)
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

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_BAM: &str = "tests/mt.sorted.bam";
    const CHR_M: &str = "chrM";

    #[test]
    fn test_sequential_count_matches_noodles() -> Result<(), Box<dyn std::error::Error>> {
        // Verify parallel BAI-interval count equals sequential noodles reader count.
        let noodles_seq = crate::count_from_standard_bam_reader(TEST_BAM, 1)? as usize;
        let bai_path = format!("{}.bai", TEST_BAM);
        let linear_indexes = get_linear_indexes(&bai_path)?;
        let intervals = get_linear_intervals(&linear_indexes)?;
        let all = crate::bam_parser::get_entire_bam_intervals(TEST_BAM, &intervals)?;
        let chunks = crate::bam_parser::merge_intervals(&all, rayon::current_num_threads().max(1));
        let par_count: usize = chunks
            .into_par_iter()
            .map(|(s, e)| crate::bam_parser::read_bam_by_interval(TEST_BAM, s, e)
                .map(|mut r| r.records().count()))
            .collect::<io::Result<Vec<_>>>()?
            .into_iter().sum();
        assert_eq!(par_count, noodles_seq,
            "parallel interval count {par_count} != sequential noodles {noodles_seq}");
        Ok(())
    }

    #[test]
    fn test_fetch_chrm_count_matches_sequential() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_sam::alignment::Record as _;

        let header = crate::get_bam_header(TEST_BAM)?;
        let bai_path = format!("{}.bai", TEST_BAM);
        let fetched = fetch_contig_or_region(TEST_BAM, &bai_path, &header, CHR_M, None, None)?;
        assert!(!fetched.is_empty(), "expected records on {CHR_M}");

        // Sequential reader filtered to chrM as ground truth.
        let (chrm_id, _, _) = header.reference_sequences()
            .get_full(CHR_M.as_bytes())
            .ok_or("chrM not found in header")?;
        let file = File::open(TEST_BAM)?;
        let mut reader = noodles_bam::io::Reader::new(file);
        let _ = reader.read_header()?;
        let seq_count = reader.records()
            .filter_map(|r| r.ok())
            .filter(|rec| {
                rec.reference_sequence_id().transpose().ok().flatten()
                    .map(|id| id == chrm_id)
                    .unwrap_or(false)
            })
            .count();

        assert_eq!(
            fetched.len(), seq_count,
            "fetch_contig_or_region count {} != sequential filter count {}",
            fetched.len(), seq_count
        );
        Ok(())
    }

    #[test]
    fn test_mapped_record_fields_valid() -> Result<(), Box<dyn std::error::Error>> {
        let header = crate::get_bam_header(TEST_BAM)?;
        let bai_path = format!("{}.bai", TEST_BAM);
        let recs = fetch_contig_or_region(TEST_BAM, &bai_path, &header, CHR_M, None, None)?;
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
                assert!(q.iter().all(|&v| v <= 93), "Phred score > 93");
            }
        }
        Ok(())
    }

    #[test]
    fn test_fetch_region_subset_and_bounds() -> Result<(), Box<dyn std::error::Error>> {
        let header = crate::get_bam_header(TEST_BAM)?;
        let bai_path = format!("{}.bai", TEST_BAM);
        let all    = fetch_contig_or_region(TEST_BAM, &bai_path, &header, CHR_M, None, None)?;
        let stop   = 2_000i64;
        let region = fetch_contig_or_region(TEST_BAM, &bai_path, &header, CHR_M, Some(0), Some(stop))?;

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

    #[test]
    fn test_flag_getters_consistent() -> Result<(), Box<dyn std::error::Error>> {
        let header = crate::get_bam_header(TEST_BAM)?;
        let bai_path = format!("{}.bai", TEST_BAM);
        let recs = fetch_contig_or_region(TEST_BAM, &bai_path, &header, CHR_M, None, None)?;
        assert!(!recs.is_empty());
        for d in &recs {
            let f = d.flag;
            let is_rev = f & 0x010 != 0;
            let is_fwd = f & 0x010 == 0;
            assert_ne!(is_rev, is_fwd, "flag 0x{f:04x}: reverse and forward must differ");

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
                tags:                 vec![],
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

    #[test]
    fn test_header_chrm_length() -> Result<(), Box<dyn std::error::Error>> {
        let header = crate::get_bam_header(TEST_BAM)?;
        assert!(header.reference_sequences().len() > 0);
        let (_, _, chrm_seq) = header.reference_sequences()
            .get_full(CHR_M.as_bytes())
            .ok_or("chrM not found in header")?;
        assert_eq!(chrm_seq.length().get(), 16569, "chrM length should be 16569 bp");
        Ok(())
    }

    // Noodles parallel fetch must return the same total count as the sequential reader.
    #[test]
    fn test_noodles_fetch_count_matches_htslib() -> Result<(), Box<dyn std::error::Error>> {
        use crate::bai_parser::{get_linear_indexes, get_linear_intervals};
        use crate::bam_parser::{get_entire_bam_intervals, merge_intervals, read_bam_by_interval};

        let bai_path = format!("{}.bai", TEST_BAM);
        let header = crate::get_bam_header(TEST_BAM)?;

        let linear_indexes = get_linear_indexes(&bai_path)?;
        let intervals = get_linear_intervals(&linear_indexes)?;
        let all = get_entire_bam_intervals(TEST_BAM, &intervals)?;
        let chunks = merge_intervals(&all, rayon::current_num_threads().max(1));

        let noodles_count: usize = chunks
            .into_par_iter()
            .map(|(start, end)| -> io::Result<usize> {
                let mut reader = read_bam_by_interval(TEST_BAM, start, end)?;
                let mut n = 0usize;
                for result in reader.records() {
                    RecordData::from_noodles(&result?, &header)?;
                    n += 1;
                }
                Ok(n)
            })
            .collect::<io::Result<Vec<usize>>>()
            .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)?
            .into_iter()
            .sum();

        let bai_path2 = format!("{}.bai", TEST_BAM);
        let chrm_count = fetch_contig_or_region(TEST_BAM, &bai_path2, &header, CHR_M, None, None)?.len();

        // noodles_count covers all chromosomes; chrm_count is chrM only.
        assert!(
            noodles_count >= chrm_count,
            "noodles total {noodles_count} should be >= chrM count {chrm_count}"
        );

        // Also verify against noodles standard reader.
        let standard = crate::count_from_standard_bam_reader(TEST_BAM, 1)? as usize;
        assert_eq!(
            noodles_count, standard,
            "noodles parallel fetch count {noodles_count} != standard reader {standard}"
        );

        Ok(())
    }

    // from_noodles must populate basic fields correctly for mapped reads.
    #[test]
    fn test_from_noodles_fields_valid() -> Result<(), Box<dyn std::error::Error>> {
        use crate::bai_parser::{get_linear_indexes, get_linear_intervals};
        use crate::bam_parser::{get_entire_bam_intervals, merge_intervals, read_bam_by_interval};

        let bai_path = format!("{}.bai", TEST_BAM);
        let header = crate::get_bam_header(TEST_BAM)?;

        let linear_indexes = get_linear_indexes(&bai_path)?;
        let intervals = get_linear_intervals(&linear_indexes)?;
        let all = get_entire_bam_intervals(TEST_BAM, &intervals)?;
        let chunks = merge_intervals(&all, 1);

        let mut checked = 0u32;
        'outer: for (start, end) in chunks {
            let mut reader = read_bam_by_interval(TEST_BAM, start, end)?;
            for result in reader.records() {
                let rec = result?;
                let d = RecordData::from_noodles(&rec, &header)?;
                if d.flag & 0x004 != 0 { continue; } // skip unmapped
                assert!(d.reference_id >= 0, "reference_id < 0 for mapped read");
                assert!(d.reference_start >= 0, "reference_start < 0 for mapped read");
                assert!(!d.query_sequence.is_empty(), "empty sequence");
                assert!(!d.cigartuples.is_empty(), "empty cigar for mapped read");
                checked += 1;
                if checked >= 100 { break 'outer; }
            }
        }
        assert!(checked > 0, "no mapped records checked");
        Ok(())
    }

    #[test]
    fn test_idxstats_basic() -> Result<(), Box<dyn std::error::Error>> {
        use noodles::csi::BinningIndex;
        use noodles::csi::binning_index::ReferenceSequence as _;

        let bai_path = format!("{}.bai", TEST_BAM);
        let index = noodles_bam::bai::fs::read(&bai_path)?;
        let header = crate::get_bam_header(TEST_BAM)?;

        // Collect idxstats directly via the index API.
        let stats: Vec<(String, u64, u64, u64)> = header.reference_sequences().keys()
            .map(|k| String::from_utf8_lossy(k).into_owned())
            .zip(header.reference_sequences().values().map(|v| v.length().get() as u64))
            .zip(index.reference_sequences())
            .map(|((name, len), ref_seq)| {
                let (n_mapped, n_unmapped) = ref_seq.metadata()
                    .map(|m| (m.mapped_record_count(), m.unmapped_record_count()))
                    .unwrap_or((0, 0));
                (name, len, n_mapped, n_unmapped)
            })
            .collect();

        // chrM must be present with non-zero mapped count.
        let chrm = stats.iter().find(|(name, _, _, _)| name == CHR_M)
            .ok_or("chrM not found in idxstats")?;
        assert!(chrm.2 > 0, "chrM n_mapped should be > 0");
        assert_eq!(chrm.1, 16569, "chrM length should be 16569");

        // BAI invariant: n_mapped + n_unmapped_per_ref + unplaced = total records.
        let total_mapped:   u64 = stats.iter().map(|(_, _, m, _)| m).sum();
        let total_unmapped: u64 = stats.iter().map(|(_, _, _, u)| u).sum();
        let unplaced = index.unplaced_unmapped_record_count().unwrap_or(0);
        let idxstats_total = total_mapped + total_unmapped + unplaced;
        let sequential = crate::count_from_standard_bam_reader(TEST_BAM, 1)?;
        assert_eq!(
            idxstats_total, sequential,
            "idxstats total (mapped+unmapped+unplaced={idxstats_total}) != sequential count {sequential}"
        );

        Ok(())
    }
}
