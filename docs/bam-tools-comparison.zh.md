# BAM Parallel Reader Comparison: bamstorm vs RabbitBAM vs quickbam

## Overview

| | bamstorm | RabbitBAM | quickbam |
|---|---|---|---|
| Language | Rust | C++ | C++ |
| Parallelism | rayon (work-stealing) | pthreads (manual) | Intel TBB |
| Decompressor | libdeflate (via noodles-bgzf) | libdeflate | libdeflate |
| Index required | Yes (BAI) | Yes (BAI) | Optional |
| License | — | — | — |

---

## Threading Model

### bamstorm

```
rayon thread pool
  ├── worker 0: seek → read interval → decompress → count
  ├── worker 1: seek → read interval → decompress → count
  └── worker N: seek → read interval → decompress → count
```

- BAI linear index intervals merged into `N=thread_count` chunks（`merge_intervals`）
- 每个 rayon worker 独立完成完整链路：I/O、解压、record 扫描串行执行
- I/O 与解压**无重叠**（串行链路）
- 高线程数下扩展性好，无跨线程 bottleneck

### RabbitBAM

```
read_thread (1)  →  [BamRead queue]  →  compress_thread × N  →  [BamCompress queue]  →  assign_thread (1)
     I/O                                      解压                                           解析/计数
```

- 严格三阶段流水线，环形缓冲 + 互斥锁 + 忙轮询（`sleep_for 1-5ms`）
- I/O 与解压**并行重叠**，是低线程数的主要性能优势
- 单线程 I/O + 单线程 parse 是高并发瓶颈（128 线程时退化至 41s，不如 bamstorm 的 3.8s）
- `pthread_setaffinity_np` 绑核，减少上下文切换

### quickbam

```
TBB parallel_for (regions)
  └── worker K:
        slice(coffset_start, coffset_end)  ← pread 整段压缩数据到内存
        bgzf_inflate_range_p(buffer)       ← TBB 并行解压 buffer 内各 block
        nfo_iterator 扫描解压数据           ← 仅读 4-byte block_size，跳过 record body
```

- 按 BAI interval 切 region，region 级并行（TBB）
- 每个 region 先整体读入内存再解压，I/O 与解压**无重叠**
- 额外支持无索引模式：将文件切 10MB chunk，启发式扫描 BGZF magic + `bam_is_valid` 定位 record 起点

---

## I/O Strategy

| | bamstorm | RabbitBAM | quickbam |
|---|---|---|---|
| 读取方式 | `File::open` + `read_exact` per interval | 专用 I/O 线程逐 block `fread` | `pread` 整段到内存 buffer |
| I/O 并行度 | 多线程各自独立读 | 单线程顺序读 | 多线程各自 pread（region 粒度）|
| I/O/解压重叠 | ❌ | ✅（流水线） | ❌ |
| mmap | ❌ | ❌ | ✅（`mfile_t`，zero-copy） |
| 内存峰值 | 一次读入所有 interval | 环形缓冲（固定大小） | 一次读入一个 region |

---

## Record Counting Hot Path

### bamstorm（当前）

```rust
loop {
    bgzf.read_exact(&mut block_size_buf)?;         // 读 4 字节
    let block_size = u32::from_le_bytes(...);
    record_data.resize(block_size, 0);
    bgzf.read_exact(&mut record_data)?;            // 读整个 record body（仅为跳过）
    n += 1;
}
```

**问题**：record body 全部解压后还需读入内存，浪费内存带宽。

### quickbam

```cpp
// nfo_iterator.h
next_ptr = cur_ptr + sizeof(uint32_t) + *((uint32_t*)cur_ptr);
// 读 4 字节，指针直接跳到下一条，零拷贝，零分支
```

**优势**：直接在解压后的内存上做指针运算，不触碰 record body。

---

## Decompression

三者均使用 libdeflate，核心解压性能相当。差异在于调用方式：

| | bamstorm | RabbitBAM | quickbam |
|---|---|---|---|
| 调用路径 | noodles-bgzf → libdeflate | 直接调用 libdeflate C API | 直接调用 libdeflate C API |
| 解压粒度 | 每个 worker 内串行按 block 解压 | 多线程 worker 并行按 block | TBB 并行按 block（region 内）|
| buffer 策略 | `Cursor<Vec<u8>>` per interval | 预分配环形缓冲池 | pread 到堆上 `unique_ptr<uint8_t[]>` |

---

## Index Usage

| | bamstorm | RabbitBAM | quickbam |
|---|---|---|---|
| 索引格式 | BAI | BAI | BAI（可选）|
| 无索引支持 | ❌ | ❌ | ✅（启发式 BGZF 扫描）|
| 分块策略 | linear index intervals → merge 成 N chunks | linear index intervals | linear index intervals → regions |
| interval 合并 | ✅（`merge_intervals`，按压缩字节均分） | ❌ | ❌ |

---

## Benchmark Results（899M records，128 cores）

| threads | bamstorm | RabbitBAM | samtools |
|---------|----------|-----------|---------|
| 1 | 270s / 181 MB/s | 231s / 212 MB/s | 253s / 193 MB/s |
| 2 | 128s / 383 MB/s | 117s / 419 MB/s | — |
| 4 | 68s / 721 MB/s | 59s / 827 MB/s | — |
| 8 | 36s / 1375 MB/s | 31s / 1598 MB/s | — |
| 128 | **3.8s / 12872 MB/s** | 41s / 1189 MB/s | 61s / 803 MB/s |

- 低线程（1-8）：RabbitBAM 领先约 15%，来源于流水线 I/O 重叠 + record body 跳过
- 高线程（128）：bamstorm 领先 10x，RabbitBAM 单线程 I/O 和 parse 成为瓶颈

---

## Pros & Cons

### bamstorm

**优点**
- 高线程扩展性极好，无单点 bottleneck
- 代码简单，rayon 自动负载均衡
- Rust 内存安全，无 data race 风险
- `merge_intervals` 减少 seek 次数，低线程下有额外收益

**缺点**
- record body 全读（count 时浪费内存带宽）
- I/O 与解压无流水线重叠
- 必须有 BAI 索引

---

### RabbitBAM

**优点**
- 三阶段流水线让 I/O 与解压重叠，低线程吞吐好
- 预分配缓冲池，内存分配稳定

**缺点**
- 单线程 I/O + 单线程 parse 限制高并发扩展性（128 线程退化）
- 忙轮询（`sleep_for`）浪费 CPU
- 手动线程管理复杂，代码维护成本高

---

### quickbam

**优点**
- record 计数热路径最优（指针跳转，零拷贝）
- 支持无索引并行读取（唯一）
- mmap zero-copy，减少一次内存拷贝
- TBB work-stealing，负载均衡好

**缺点**
- 先整段 pread 进内存再解压，region 较大时内存峰值高
- I/O 与解压无流水线重叠
- 依赖 Intel TBB，非标准依赖

---

## Key Takeaways

1. **流水线 vs 简单并行**：RabbitBAM 的流水线在低线程时有优势，但单线程 I/O/parse 成为高并发天花板；bamstorm 和 quickbam 的"每 worker 自治"模型扩展性更好。

2. **record body 跳过**：quickbam 的 `nfo_iterator` 是低线程数下最有效的单点优化，bamstorm 引入后预计可消除与 RabbitBAM 的 15% 差距。

3. **无索引支持**：quickbam 独有，适合处理未建索引的 BAM 文件，是功能层面的差异化。

4. **高并发场景**：bamstorm 在 128 线程下以 12872 MB/s 遥遥领先，适合大规模并行计算环境。
