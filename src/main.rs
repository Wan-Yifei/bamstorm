use bamstrom::{
    bai_parser::get_linear_indexes, bai_parser::get_linear_intervals, get_bam_header,
    process_records_through_intervals, update_bam_record,
};
use noodles::bam as noodles_bam;
use std::error::Error;
use std::{collections::HashMap, fs::File};

fn main() -> Result<(), Box<dyn Error>> {
    let output_path = "output.bam";
    let input_path = "tests/chr1.bam";
    let input_bai_path = "tests/chr1.bam.bai";

    let header = get_bam_header(input_path)?;
    let linear_indexes_all = get_linear_indexes(input_bai_path)?;
    let intervals = get_linear_intervals(&linear_indexes_all)?;

    let mut updated_fields = HashMap::new();
    updated_fields.insert("name", "NEW_NAME");

    let mut writer = noodles_bam::io::Writer::new(File::create(output_path)?);
    writer.write_header(&header)?;

    process_records_through_intervals(
        input_path,
        &intervals,
        &header,
        writer,
        |record| {
            let updated =
                update_bam_record(&header, record, &updated_fields, "SOLEXA-1GA-2_2_FC20EMB")?;
            Ok(Some(updated))
        },
    )?;

    Ok(())
}
