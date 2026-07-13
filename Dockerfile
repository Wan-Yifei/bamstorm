# Stage 1: build bamstorm binaries + Python extension wheel
FROM rust:slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    pkg-config \
    libhts-dev \
    libclang-dev \
    python3 \
    python3-dev \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install maturin --break-system-packages

WORKDIR /build
COPY Cargo.toml Cargo.lock pyproject.toml ./
COPY src ./src
COPY python ./python

RUN cargo build --release --bin bench_count
RUN maturin build --release --out /tmp/wheels

# Stage 2: build RabbitBAM and its dependencies (htslib + libdeflate) from source
FROM debian:bookworm-slim AS rabbitbam-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    make \
    g++ \
    autoconf \
    automake \
    libtool \
    pkg-config \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# libdeflate >= 1.12 (RabbitBAM dependency)
RUN git clone --depth 1 --branch v1.22 https://github.com/ebiggers/libdeflate.git /opt/libdeflate \
    && cmake -S /opt/libdeflate -B /opt/libdeflate/build -DCMAKE_INSTALL_PREFIX=/usr/local \
    && cmake --build /opt/libdeflate/build -j$(nproc) \
    && cmake --install /opt/libdeflate/build

# htslib >= 1.15 (RabbitBAM dependency)
RUN git clone --depth 1 --branch 1.21 https://github.com/samtools/htslib.git /opt/htslib \
    && cd /opt/htslib \
    && git submodule update --init --recursive \
    && autoreconf -i \
    && ./configure --prefix=/usr/local --with-libdeflate \
    && make -j$(nproc) \
    && make install

# RabbitBAM — patch config.h to fix narrowing conversion error under GCC 12
RUN git clone https://github.com/RabbitBio/RabbitBAM.git /opt/RabbitBAM \
    && sed -i 's/{-1, 0, -1, 1, 2, -1, -1, 3}/{(int8)-1, 0, (int8)-1, 1, 2, (int8)-1, (int8)-1, 3}/' /opt/RabbitBAM/config.h

WORKDIR /opt/RabbitBAM
SHELL ["/bin/bash", "-c"]

RUN bash configure.sh /usr/local /usr/local \
    && source env.sh \
    && make clean \
    && make -j$(nproc)

# Stage 3: runtime image with samtools + python3-pysam + RabbitBAM
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    python3 \
    python3-pip \
    python3-pysam \
    fio \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# bamstorm CLI counter
COPY --from=builder /build/target/release/bench_count /app/bench_count

# bamstorm Python extension wheel
COPY --from=builder /tmp/wheels /tmp/wheels
RUN pip3 install /tmp/wheels/bamstorm-*.whl --break-system-packages \
    && rm -rf /tmp/wheels

# RabbitBAM binary + all build artifacts (binary links against tools.o and .so files)
COPY --from=rabbitbam-builder /opt/RabbitBAM /opt/RabbitBAM

# htslib / libdeflate shared libs needed at runtime
COPY --from=rabbitbam-builder /usr/local/lib/libhts.so* /usr/local/lib/
COPY --from=rabbitbam-builder /usr/local/lib/libdeflate.so* /usr/local/lib/

RUN ldconfig

ENV LD_LIBRARY_PATH=/opt/RabbitBAM:/usr/local/lib

# Benchmark scripts and config
COPY bench/bench.py /app/bench.py
COPY bench/bench_coverage.py /app/bench_coverage.py
COPY bench/bench.toml /app/bench.toml

# Test data is mounted at runtime; nothing to COPY here.
# Usage:
#   docker run --rm -v /path/to/data:/data bamstorm-bench \
#       python3 /app/bench.py /data/sample.bam /data/sample.bam.bai
CMD ["python3", "/app/bench.py", "--help"]
