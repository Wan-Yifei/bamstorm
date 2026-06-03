# BAM Parallel Reader Comparison: bamstrom vs RabbitBAM vs quickbam

## Overview

| | bamstrom | RabbitBAM | quickbam |
|---|---|---|---|
| Language | Rust | C++ | C++ |
| Parallelism | rayon (work-stealing) | pthreads (manual pipeline) | Intel TBB |
| Decompressor | libdeflate (via noodles-bgzf) | libdeflate | libdeflate |
| BAI index required | Yes | Yes | Optional |

---

## Threading Model

### bamstrom

```
rayon thread pool
  ├── worker 0: seek → read interval → decompress → count
  ├── worker 1: seek → read interval → decompress → count
  └── worker N: seek → read interval → decompress → count
```

BAI linear index intervals are merged into `N = thread_count` equal-sized chunks
(`merge_intervals`). Each rayon worker owns its full pipeline — I/O, decompression,
and record scanning are serialized within a worker. No inter-thread synchronization
is needed, which gives excellent scaling at high thread counts.

### RabbitBAM

```
read_thread (×1) → [BamRead queue] → compress_thread (×N) → [BamCompress queue] → assign_thread (×1)
      I/O                                  decompress                                   parse / count
```

A strict three-stage pipeline backed by ring buffers, mutexes, and busy-polling
(`sleep_for 1–5 ms`). I/O and decompression overlap, which helps at low thread
counts. However, the single-threaded I/O and single-threaded parse stages become
bottlenecks at high concurrency (degraded to 41 s at 128 threads vs bamstrom's 3.8 s).
CPU affinity is pinned per thread via `pthread_setaffinity_np`.

### quickbam

```
tbb::parallel_for (over regions)
  └── worker K:
        slice(coffset_start, coffset_end)   ← pread entire compressed region into memory
        bgzf_inflate_range_p(buffer, ...)   ← TBB parallel decompress blocks in buffer
        nfo_iterator scan                   ← read 4-byte block_size, skip record body
```

Regions are derived from BAI intervals. Each TBB worker reads its own region with a
single `pread` call, decompresses all BGZF blocks inside in parallel (TBB), then
scans the decompressed buffer via pointer arithmetic. I/O and decompression do not
overlap within a region. Uniquely supports index-free parallel reading by splitting
the file into 10 MB chunks and heuristically scanning for BGZF magic bytes.

---

## I/O Strategy

| | bamstrom | RabbitBAM | quickbam |
|---|---|---|---|
| Read mechanism | `File::open` + `read_exact` per interval | Dedicated I/O thread, block-by-block `fread` | `pread` entire region into heap buffer |
| I/O parallelism | Each worker reads independently | Single thread, sequential | Each worker `pread`s its own region |
| I/O / decompress overlap | No | **Yes (pipeline)** | No |
| mmap | No | No | **Yes** (`mfile_t`, zero-copy) |
| Memory peak | All intervals loaded simultaneously | Fixed-size ring buffer | One region at a time per worker |

---

## Record Counting Hot Path

### bamstrom (current)

```rust
loop {
    bgzf.read_exact(&mut block_size_buf)?;      // read 4-byte length prefix
    let block_size = u32::from_le_bytes(...);
    record_data.resize(block_size, 0);
    bgzf.read_exact(&mut record_data)?;         // read full record body — only to discard it
    n += 1;
}
```

The entire record body is decompressed and copied into memory just to be skipped,
wasting memory bandwidth proportional to the average record size.

### quickbam

```cpp
// nfo_iterator.h
next_ptr = cur_ptr + sizeof(uint32_t) + *((uint32_t*)cur_ptr);
```

A single pointer arithmetic expression advances past each record. No data is read
beyond the 4-byte length prefix. Zero copies, zero branches — the tightest possible
inner loop for BAM record counting.

---

## Decompression

All three tools use libdeflate. The differences lie in how it is invoked:

| | bamstrom | RabbitBAM | quickbam |
|---|---|---|---|
| Call path | noodles-bgzf → libdeflate C API | Direct libdeflate C API | Direct libdeflate C API |
| Granularity | Per-block, serialized within each worker | Per-block, parallel across pipeline workers | Per-block, TBB parallel within each region |
| Buffer strategy | `Cursor<Vec<u8>>` per interval | Pre-allocated ring buffer pool | Heap-allocated `unique_ptr<uint8_t[]>` per region |

---

## Index Usage

| | bamstrom | RabbitBAM | quickbam |
|---|---|---|---|
| Index format | BAI | BAI | BAI (optional) |
| Index-free reading | No | No | **Yes** (heuristic BGZF scan) |
| Chunking strategy | Linear index intervals → merged into N chunks | Linear index intervals | Linear index intervals → regions |
| Interval merging | **Yes** (`merge_intervals`, balanced by compressed bytes) | No | No |

---

## Benchmark Results (899 M records, 128-core machine)

| Threads | bamstrom | RabbitBAM | samtools |
|---------|----------|-----------|---------|
| 1 | 270 s / 181 MB/s | 231 s / 212 MB/s | 253 s / 193 MB/s |
| 2 | 128 s / 383 MB/s | 117 s / 419 MB/s | — |
| 4 | 68 s / 721 MB/s | 59 s / 827 MB/s | — |
| 8 | 36 s / 1375 MB/s | 31 s / 1598 MB/s | — |
| 128 | **3.8 s / 12872 MB/s** | 41 s / 1189 MB/s | 61 s / 803 MB/s |

- **Low thread counts (1–8):** RabbitBAM leads by ~15%, primarily due to I/O–decompress
  pipeline overlap and body-skipping in the record scan.
- **High thread counts (128):** bamstrom leads by ~10×. RabbitBAM's single-threaded I/O
  and single-threaded parse stages are fully saturated, starving the decompression workers.

---

## Pros & Cons

### bamstrom

**Strengths**
- Excellent scaling at high thread counts — no single-threaded bottleneck
- Simple architecture; rayon handles load balancing automatically
- Memory-safe by construction (Rust); no data races possible
- `merge_intervals` reduces seek overhead, improving low-thread throughput

**Weaknesses**
- Reads full record body during counting (wasted memory bandwidth)
- No I/O / decompress pipeline overlap
- Requires BAI index

---

### RabbitBAM

**Strengths**
- Pipeline overlap between I/O and decompression improves single-thread throughput
- Pre-allocated ring buffer pool keeps memory allocation deterministic

**Weaknesses**
- Single-threaded I/O and single-threaded parse cap scaling at high concurrency
- Busy-polling (`sleep_for`) wastes CPU cycles while waiting
- Manual thread management increases code complexity and maintenance cost

---

### quickbam

**Strengths**
- Tightest record-counting inner loop (pointer arithmetic, zero copies)
- Only tool supporting index-free parallel reading
- mmap zero-copy eliminates one memory copy per region
- TBB work-stealing provides automatic load balancing

**Weaknesses**
- Entire region loaded into memory before decompression begins; memory peak scales with region size
- No I/O / decompress pipeline overlap
- Dependency on Intel TBB

---

## Key Takeaways

1. **Pipeline vs. autonomous workers.** RabbitBAM's pipeline gives an edge at low
   thread counts by overlapping I/O and decompression. At high thread counts its
   single-threaded I/O and parse stages become the ceiling. bamstrom and quickbam's
   per-worker autonomous model scales linearly.

2. **Skipping record bodies.** quickbam's `nfo_iterator` is the most impactful
   single optimization for counting workloads. Applying it to bamstrom is expected
   to close most of the remaining ~15% gap against RabbitBAM at low thread counts.

3. **Index-free reading.** quickbam's heuristic BGZF scan is a functional
   differentiator — the only option when BAI is unavailable.

4. **High-concurrency workloads.** bamstrom at 128 threads (12,872 MB/s) outperforms
   both alternatives by a wide margin, making it the best choice for large-scale
   parallel compute environments.
