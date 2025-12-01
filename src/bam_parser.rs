use std::{
    fs::File,
    io::{self, BufWriter},
};
use noodles::bam as bam;
use noodles::sam::{self as sam, alignment::io::Write};

pub fn standard_bam_read() -> io::Result<()> {
    let src = "e:/project/bamstrom/tests/mt.sorted.bam"; 

    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    // let stdout = io::stdout().lock();
    // let mut writer = sam::io::Writer::new(BufWriter::new(stdout));
    let file = File::create("test.sam")?;
    let mut writer = sam::io::Writer::new(file);
    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[ignore]
    fn test_standard_bam_read() -> Result<(), Box<dyn std::error::Error>> {
        // Run bam_read and propagate any errors
        standard_bam_read()?;  

        Ok(())
    }
}
