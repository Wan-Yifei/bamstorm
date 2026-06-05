/// De-identify a BAM file by replacing a PHI pattern in read names.
///
/// Usage:
///   cargo run --example deidentify -- \
///       <input.bam> <input.bai> <output.bam> <phi_pattern> <replacement>
///
/// Example:
///   cargo run --example deidentify -- \
///       patient.bam patient.bam.bai anon.bam "PATIENT-001" "SUBJECT"
use bamstorm::{
    bai_parser::{get_linear_indexes, get_linear_intervals},
    get_bam_header, process_records_through_intervals, update_bam_record,
};
use noodles::bam as noodles_bam;
use std::{collections::HashMap, error::Error, fs::File};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 6 {
        eprintln!(
            "Usage: {} <input.bam> <input.bai> <output.bam> <phi_pattern> <replacement>",
            args[0]
        );
        std::process::exit(1);
    }
    let (input_bam, input_bai, output_bam) = (&args[1], &args[2], &args[3]);
    let (phi_pattern, replacement) = (&args[4], &args[5]);

    let header = get_bam_header(input_bam)?;
    let intervals = get_linear_intervals(&get_linear_indexes(input_bai)?)?;

    let mut updated_fields = HashMap::new();
    updated_fields.insert("name", replacement.as_str());

    let mut writer = noodles_bam::io::Writer::new(File::create(output_bam)?);
    writer.write_header(&header)?;

    let n = std::sync::atomic::AtomicU64::new(0);
    process_records_through_intervals(
        input_bam,
        &intervals,
        &header,
        writer,
        |record| {
            n.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let updated = update_bam_record(&header, record, &updated_fields, phi_pattern)?;
            Ok(Some(updated))
        },
    )?;

    println!(
        "Processed {} records. Output written to {output_bam}",
        n.load(std::sync::atomic::Ordering::Relaxed)
    );
    Ok(())
}
