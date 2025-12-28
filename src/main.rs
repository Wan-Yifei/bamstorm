use bamstrom::{get_bam_header, update_bam_record, bai_parser::get_linear_indexes, bai_parser::get_linear_intervals, process_records_through_intervals};
use std::sync::{Arc, Mutex};
use std::{fs::File, collections::HashMap};
use std::error::Error;
use noodles::bam as noodles_bam;

fn main() -> Result<(), Box<dyn Error>> {
    let output_path = "output.bam";
    let mut writer = noodles_bam::io::Writer::new(File::create(output_path).unwrap());

    let input_path = "tests/chr1.bam";
    let input_bai_path = "tests/chr1.bam.bai";
    // Set Multithread Reader
    let _header = get_bam_header(input_path)?;

    // Get linear intervals from BAI

    let linear_indexes_all = get_linear_indexes(input_bai_path)?;
    let intervals = get_linear_intervals(&linear_indexes_all)?;

    let mut updated_fields = HashMap::new();
    updated_fields.insert("name", "NEW_NAME");

    // Output BAM
    writer.write_header(&_header)?;

    // Process BAM records
    let output_writer = Arc::new(Mutex::new(writer));
    process_records_through_intervals(input_path, &intervals, 4, output_writer, |record, output_writer| {
        // For demonstration, replace read name old_name with "new_name"
        let updated_record = update_bam_record(&_header, &record, &updated_fields, "SOLEXA-1GA-2_2_FC20EMB", output_writer)?;
        Ok(updated_record)
    })?;


    Ok(())
}

