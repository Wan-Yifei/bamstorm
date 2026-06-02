use bamstrom::bai_parser::{get_linear_indexes, get_linear_intervals};
use bamstrom::bam_parser::{count_records_in_virtual_range, get_entire_bam_intervals};
use rayon::prelude::*;
use std::io;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();

    let mut bam_path: Option<&str> = None;
    let mut bai_path: Option<&str> = None;
    let mut threads: usize = 0;
    let mut verbose = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--threads" | "-t" => {
                i += 1;
                threads = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(0);
            }
            "--verbose" | "-v" => verbose = true,
            arg if bam_path.is_none() => bam_path = Some(arg),
            arg if bai_path.is_none() => bai_path = Some(arg),
            _ => {}
        }
        i += 1;
    }

    let (bam_path, bai_path) = match (bam_path, bai_path) {
        (Some(b), Some(i)) => (b, i),
        _ => {
            eprintln!("Usage: bench_count [--threads N] [--verbose] <bam_path> <bai_path>");
            std::process::exit(1);
        }
    };

    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap_or(());
    }

    let linear_indexes = get_linear_indexes(bai_path)?;
    let intervals = get_linear_intervals(&linear_indexes)?;
    let all_intervals = get_entire_bam_intervals(bam_path, &intervals)?;

    if verbose {
        // Check for gap before the first interval (reads between header and first index entry)
        eprintln!(
            "[verbose] {} intervals  ({}  BAI-derived + 1 tail)",
            all_intervals.len(),
            all_intervals.len() - 1,
        );
        eprintln!(
            "[verbose] first interval start: compressed={:#x} uncompressed={}",
            all_intervals[0].0.compressed(),
            all_intervals[0].0.uncompressed(),
        );
        eprintln!(
            "[verbose] last  interval end  : compressed={:#x} uncompressed={}",
            all_intervals.last().unwrap().1.compressed(),
            all_intervals.last().unwrap().1.uncompressed(),
        );

        // Per-interval counts (sequential so output is ordered)
        let mut total = 0u64;
        for (idx, &(start, end)) in all_intervals.iter().enumerate() {
            let n = count_records_in_virtual_range(bam_path, start, end)?;
            total += n;
            if n == 0 {
                eprintln!(
                    "[verbose] interval {:4}: 0 records  [{:#x}:{} → {:#x}:{}]",
                    idx,
                    start.compressed(), start.uncompressed(),
                    end.compressed(),   end.uncompressed(),
                );
            }
        }
        eprintln!("[verbose] total = {total}");
        println!("{total}");
        return Ok(());
    }

    let count: u64 = all_intervals
        .into_par_iter()
        .map(|(start, end)| count_records_in_virtual_range(bam_path, start, end))
        .sum::<io::Result<u64>>()?;

    println!("{count}");
    Ok(())
}
