FROM rust:1.92

# Install htop and vim
RUN apt-get update && apt-get install -y \
    htop \
    vim \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /bamstorm

# Copy entire project
COPY ./src src
COPY Cargo.lock Cargo.lock
COPY Cargo.toml Cargo.toml

# Default to interactive bash
CMD ["/bin/bash"]