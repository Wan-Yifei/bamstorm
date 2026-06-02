use bamstrom::bai_parser::{get_linear_indexes, get_linear_intervals};
use bamstrom::bam_parser::{get_entire_bam_intervals, read_bam_by_interval};
use rayon::prelude::*;
use std::io;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();

    let mut bam_path: Option<&str> = None;
    let mut bai_path: Option<&str> = None;
    let mut threads: usize = 0; // 0 = rayon default (all cores)

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--threads" | "-t" => {
                i += 1;
                threads = args
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
            }
            arg if bam_path.is_none() => bam_path = Some(arg),
            arg if bai_path.is_none() => bai_path = Some(arg),
            _ => {}
        }
        i += 1;
    }

    let (bam_path, bai_path) = match (bam_path, bai_path) {
        (Some(b), Some(i)) => (b, i),
        _ => {
            eprintln!("Usage: bench_count [--threads N] <bam_path> <bai_path>");
            std::process::exit(1);
        }
    };

    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap_or(());
    }

    let intervals = get_linear_intervals(&get_linear_indexes(bai_path)?)?;
    let all_intervals = get_entire_bam_intervals(bam_path, &intervals)?;

    let count: u64 = all_intervals
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
        .sum::<io::Result<u64>>()?;

    println!("{count}");
    Ok(())
}
