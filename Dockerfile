FROM rust:1.92

WORKDIR /bamstorm

# Copy entire project
COPY ./src src 
COPY Cargo.lock Cargo.lock
COPY Cargo.toml Cargo.toml

