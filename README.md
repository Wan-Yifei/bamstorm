# bamstrom

A high-performance parallel BAM reader for large-scale genomics workloads, written in Rust with Python bindings via PyO3.

## How it works

Standard tools (samtools, pysam/htslib) read BAM files with a single IO thread and decompress BGZF blocks in parallel. **bamstrom** takes a different approach: it uses the BAI linear index to split the file into independent byte-range intervals, then reads and decompresses all intervals simultaneously with rayon.

```
htslib:  single fd → sequential read → parallel BGZF decompress
bamstrom: BAI intervals → N parallel fds → parallel read + decompress
```

## Installation

**Python (recommended)**

```bash
pip install bamstrom
```

Wheels are built against the stable ABI (`abi3`) and work on Python 3.8+.

**From source (requires Rust)**

```bash
git clone https://github.com/Wan-Yifei/bamstrom
cd bamstrom
pip install maturin
maturin develop --release --features python
```

## Python usage

```python
import bamstrom

# Fast parallel record count — stays entirely in Rust, no Python overhead
n = bamstrom.count("sample.bam", "sample.bam.bai")
print(f"{n} records")

# Iterate over records
with bamstrom.AlignmentFile("sample.bam", "sample.bam.bai") as af:
    for read in af:
        if read.is_unmapped:
            continue
        print(read.query_name, read.reference_start, read.cigarstring)

# .count() method on AlignmentFile
with bamstrom.AlignmentFile("sample.bam", "sample.bam.bai") as af:
    print(af.count())
```

### `BamRecord` attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `query_name` | `str \| None` | Read name |
| `flag` | `int` | SAM flag |
| `reference_id` | `int \| None` | Reference sequence index |
| `reference_start` | `int \| None` | 0-based alignment start |
| `mapping_quality` | `int \| None` | MAPQ |
| `cigarstring` | `str` | CIGAR string (e.g. `"101M"`) |
| `query_sequence` | `str` | Nucleotide sequence |
| `template_length` | `int` | TLEN |
| `is_paired` | `bool` | Flag 0x1 |
| `is_proper_pair` | `bool` | Flag 0x2 |
| `is_unmapped` | `bool` | Flag 0x4 |
| `is_mate_unmapped` | `bool` | Flag 0x8 |
| `is_reverse` | `bool` | Flag 0x10 |
| `is_secondary` | `bool` | Flag 0x100 |
| `is_qcfail` | `bool` | Flag 0x200 |
| `is_duplicate` | `bool` | Flag 0x400 |
| `is_supplementary` | `bool` | Flag 0x800 |

## Rust usage

Add to `Cargo.toml`:

```toml
[dependencies]
bamstrom = { git = "https://github.com/Wan-Yifei/bamstrom" }
```

```rust
use bamstrom::{bai_parser::{get_linear_indexes, get_linear_intervals}, count_all_records};

let indexes = get_linear_indexes("sample.bam.bai")?;
let intervals = get_linear_intervals(&indexes)?;
let total = count_all_records("sample.bam", &intervals)?;
```

## Benchmark

Run the Docker-based benchmark comparing bamstrom against samtools and pysam:

```bash
docker build -t bamstrom-bench .
docker run --rm -v /path/to/data:/data bamstrom-bench \
    python3 /app/bench.py /data/sample.bam /data/sample.bam.bai
```

Example output:

```
  Tool                          threads    elapsed    throughput  records
  ---------------------------------------------------------------------------
  bamstrom                      threads=1    3.201s    312.4 MB/s  records=45000000
  bamstrom                      threads=4    0.891s   1121.5 MB/s  records=45000000
  bamstrom                      threads=8    0.512s   1952.0 MB/s  records=45000000
  samtools view -c              threads=1    9.847s    101.5 MB/s  records=45000000
  samtools view -c              threads=8    3.201s    312.4 MB/s  records=45000000
  pysam fetch(until_eof)        threads=1   18.234s     54.8 MB/s  records=45000000
```

## Requirements

- BAM file must be coordinate-sorted and indexed (`.bai`)
- Python ≥ 3.8 (for Python bindings)
- Rust ≥ 1.85 (for building from source)
