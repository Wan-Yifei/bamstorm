use crate::{
    bai_parser::{get_linear_indexes, get_linear_intervals},
    bam_parser::{get_entire_bam_intervals, read_bam_by_interval},
    get_bam_header,
};
use bstr::ByteSlice;
use noodles::sam::{self as noodles_sam, alignment::{RecordBuf, io::Write as _}};
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use rayon::prelude::*;
use std::io;

// Accept any error type — bai functions return Box<dyn Error> or String.
fn to_py_err<E: ToString>(e: E) -> PyErr {
    PyIOError::new_err(e.to_string())
}

// Serialize one RecordBuf to SAM text then pluck the CIGAR (col 6) and
// SEQ (col 10) fields.  Slower than direct field access but avoids having
// to depend on noodles' non-Display Cigar/Sequence internals.
fn extract_cigar_seq(header: &noodles_sam::Header, rec: &RecordBuf) -> (String, String) {
    let mut buf: Vec<u8> = Vec::new();
    if noodles_sam::io::Writer::new(&mut buf)
        .write_alignment_record(header, rec)
        .is_ok()
    {
        let s = String::from_utf8_lossy(&buf);
        let cols: Vec<&str> = s.splitn(11, '\t').collect();
        return (
            cols.get(5).copied().unwrap_or("*").trim_end_matches('\n').to_string(),
            cols.get(9).copied().unwrap_or("*").trim_end_matches('\n').to_string(),
        );
    }
    ("*".into(), "*".into())
}

// ── count ─────────────────────────────────────────────────────────────────────

/// Count all records; stays entirely in Rust — no Python object overhead.
#[pyfunction]
pub fn count(bam_path: &str, bai_path: &str) -> PyResult<u64> {
    let intervals = get_linear_intervals(&get_linear_indexes(bai_path).map_err(to_py_err)?)
        .map_err(to_py_err)?;
    let all_intervals =
        get_entire_bam_intervals(bam_path, &intervals).map_err(to_py_err)?;
    all_intervals
        .into_par_iter()
        .map(|(start, end)| -> io::Result<u64> {
            let mut reader = read_bam_by_interval(bam_path, start, end)?;
            let mut n = 0u64;
            for result in reader.records() {
                result?;
                n += 1;
            }
            Ok(n)
        })
        .sum::<io::Result<u64>>()
        .map_err(to_py_err)
}

// ── internal Rust struct ──────────────────────────────────────────────────────
// Fields are extracted in Rust parallel workers (no GIL); Python objects are
// created only when the caller iterates, serialising just that step.

struct RecordData {
    query_name: Option<String>,
    flag: u16,
    reference_id: Option<i32>,
    reference_start: Option<i64>,
    mapping_quality: Option<u8>,
    cigarstring: String,
    query_sequence: String,
    template_length: i32,
}

impl RecordData {
    fn from_buf(header: &noodles_sam::Header, rec: &RecordBuf) -> Self {
        let (cigarstring, query_sequence) = extract_cigar_seq(header, rec);
        RecordData {
            query_name: rec.name().and_then(|n| n.to_str().ok()).map(String::from),
            flag: rec.flags().bits(),
            reference_id: rec.reference_sequence_id().map(|id| id as i32),
            reference_start: rec.alignment_start().map(|p| p.get() as i64 - 1),
            mapping_quality: rec.mapping_quality().map(u8::from),
            cigarstring,
            query_sequence,
            template_length: rec.template_length(),
        }
    }
}

// ── Python-visible types ──────────────────────────────────────────────────────

#[pyclass(get_all)]
pub struct BamRecord {
    pub query_name: Option<String>,
    pub flag: u16,
    pub reference_id: Option<i32>,
    pub reference_start: Option<i64>,
    pub mapping_quality: Option<u8>,
    pub cigarstring: String,
    pub query_sequence: String,
    pub template_length: i32,
}

#[pymethods]
impl BamRecord {
    fn __repr__(&self) -> String {
        let name = match &self.query_name {
            Some(n) => format!("{:?}", n),
            None => "None".to_string(),
        };
        let ref_start = match self.reference_start {
            Some(p) => p.to_string(),
            None => "None".to_string(),
        };
        format!("BamRecord(query_name={}, flag={}, ref_start={})", name, self.flag, ref_start)
    }

    #[getter] fn is_paired(&self) -> bool       { self.flag & 0x001 != 0 }
    #[getter] fn is_proper_pair(&self) -> bool   { self.flag & 0x002 != 0 }
    #[getter] fn is_unmapped(&self) -> bool      { self.flag & 0x004 != 0 }
    #[getter] fn is_mate_unmapped(&self) -> bool { self.flag & 0x008 != 0 }
    #[getter] fn is_reverse(&self) -> bool       { self.flag & 0x010 != 0 }
    #[getter] fn is_secondary(&self) -> bool     { self.flag & 0x100 != 0 }
    #[getter] fn is_qcfail(&self) -> bool        { self.flag & 0x200 != 0 }
    #[getter] fn is_duplicate(&self) -> bool     { self.flag & 0x400 != 0 }
    #[getter] fn is_supplementary(&self) -> bool { self.flag & 0x800 != 0 }
}

// ── iterator ──────────────────────────────────────────────────────────────────

#[pyclass]
pub struct RecordIterator {
    records: Vec<RecordData>,
    index: usize,
}

#[pymethods]
impl RecordIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<BamRecord> {
        let idx = slf.index;
        if idx >= slf.records.len() {
            return None;
        }
        slf.index += 1;
        let d = &slf.records[idx];
        Some(BamRecord {
            query_name: d.query_name.clone(),
            flag: d.flag,
            reference_id: d.reference_id,
            reference_start: d.reference_start,
            mapping_quality: d.mapping_quality,
            cigarstring: d.cigarstring.clone(),
            query_sequence: d.query_sequence.clone(),
            template_length: d.template_length,
        })
    }

    fn __len__(&self) -> usize {
        self.records.len().saturating_sub(self.index)
    }
}

// ── AlignmentFile ─────────────────────────────────────────────────────────────

#[pyclass]
pub struct AlignmentFile {
    bam_path: String,
    bai_path: String,
    header: noodles_sam::Header,
}

#[pymethods]
impl AlignmentFile {
    #[new]
    pub fn new(bam_path: String, bai_path: String) -> PyResult<Self> {
        let header = get_bam_header(&bam_path).map_err(to_py_err)?;
        Ok(AlignmentFile { bam_path, bai_path, header })
    }

    /// Count records without creating Python objects — fastest for pure counting.
    pub fn count(&self) -> PyResult<u64> {
        count(&self.bam_path, &self.bai_path)
    }

    /// Read all intervals in parallel; field decoding in Rust, Python objects
    /// created lazily on iteration.
    pub fn fetch(&self) -> PyResult<RecordIterator> {
        let header = &self.header;
        let bam_path = self.bam_path.as_str();
        let intervals =
            get_linear_intervals(&get_linear_indexes(&self.bai_path).map_err(to_py_err)?)
                .map_err(to_py_err)?;
        let all_intervals =
            get_entire_bam_intervals(bam_path, &intervals).map_err(to_py_err)?;

        let records: Vec<RecordData> = all_intervals
            .into_par_iter()
            .map(|(start, end)| -> io::Result<Vec<RecordData>> {
                let mut reader = read_bam_by_interval(bam_path, start, end)?;
                let mut chunk = Vec::new();
                for result in reader.records() {
                    let raw = result?;
                    let buf = RecordBuf::try_from_alignment_record(header, &raw)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                    chunk.push(RecordData::from_buf(header, &buf));
                }
                Ok(chunk)
            })
            .collect::<io::Result<Vec<_>>>()
            .map_err(to_py_err)?
            .into_iter()
            .flatten()
            .collect();

        Ok(RecordIterator { records, index: 0 })
    }

    fn __iter__(&self) -> PyResult<RecordIterator> {
        self.fetch()
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[pyo3(signature = (_exc_type=None, _exc_val=None, _exc_tb=None))]
    fn __exit__(
        &self,
        _exc_type: Option<Bound<'_, PyAny>>,
        _exc_val: Option<Bound<'_, PyAny>>,
        _exc_tb: Option<Bound<'_, PyAny>>,
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
