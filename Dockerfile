# Stage 1: build bamstrom binaries
FROM rust:slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src ./src

RUN cargo build --release --bin bench_count

# Stage 2: runtime image with samtools + python3-pysam
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    python3 \
    python3-pysam \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy compiled binary from builder
COPY --from=builder /build/target/release/bench_count /app/bench_count

# Copy benchmark script
COPY bench/bench.py /app/bench.py

# Test data is mounted at runtime; nothing to COPY here.
# Usage:
#   docker run --rm -v /path/to/data:/data bamstrom-bench \
#       python3 /app/bench.py /data/sample.bam /data/sample.bam.bai
CMD ["python3", "/app/bench.py", "--help"]
