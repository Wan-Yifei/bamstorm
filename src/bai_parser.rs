use std::fs::File;
use std::fmt::format;
use std::slice::Windows;

use noodles::bam::bai;
use noodles::bgzf::VirtualPosition;

const DEFALUT_MIN_UNCOMPRESSED_POSITION: u16 = u16::MIN;

// Get file size
fn get_file_size(file_path: &str) -> std::io::Result<u64> {
    let file = File::open(file_path)?;
    Ok(file.metadata().unwrap().len())
}

// Get all linear indexes from a BAI file
pub fn get_linear_indexes(
    bai_path: &str, bam_size: u64
) -> Result<Vec<VirtualPosition>, Box<dyn std::error::Error>> {
    // Read the BAI file
    let bai_file = bai::fs::read(bai_path)?;
    // Prepare a vector to store all linear indexes
    let mut linear_indexes: Vec<VirtualPosition> = Vec::new();

    // Iterate over all reference sequences in the BAI
    for ref_seq in bai_file.reference_sequences() {
        // Append to the master list
        linear_indexes.extend(ref_seq.index());
    }

    // Remove duplicate indexes
    linear_indexes.dedup();
    // Sort the linear indexes
    linear_indexes.sort();
    // Compare the largest linear index with the maximum compressed virtual position of BGZF file
    if linear_indexes.last().unwrap().compressed() > bam_size {
        return Err(format!("The max coffset {:?} exceeds compressed file size {:?}",
            linear_indexes.last().unwrap().compressed(), bam_size).into());
    }
    else if linear_indexes.last().unwrap().compressed() < bam_size {
        let _end_voffset = VirtualPosition::new(bam_size, DEFALUT_MIN_UNCOMPRESSED_POSITION).unwrap();
        linear_indexes.push(_end_voffset);
        println!("Added end voffset {:?} to point the file end {:?}.", _end_voffset, bam_size);
    }
    else {
        println!("All linear indexes are within valid range.");
    }

    Ok(linear_indexes)
}

pub fn convert_linear_indexes_to_coffset(linear_indexes: &[VirtualPosition]) -> Vec<u64> {
    linear_indexes.iter().map(|vp| vp.compressed()).collect()
}

pub fn get_linear_intervals(
    linear_indexes: &[VirtualPosition],
) -> Result<Vec<(VirtualPosition, VirtualPosition)>, String> {
    // Ensure there are enough indexes to form intervals
    if linear_indexes.len() < 2 {
        return Err("Not enough linear indexes to form intervals".to_string());
    }
    // Create intervals from consecutive linear indexes
    let mut intervals = Vec::new();
    for window in linear_indexes.windows(2) {
        if window[1] > window[0] {
            intervals.push((window[0], window[1]));
        } else {
            return Err(format!(
                "Linear index is not strictly increasing,
             found out-of-order or duplicate entries.: START: {:?} >  END: {:?}",
                window[0], window[1]
            ));
        }
    }
    Ok(intervals)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_bai() -> Result<(), Box<dyn std::error::Error>> {
        //let bai_path = "E:/project/bamstrom/tests/mt.sorted.bam.bai";
        // let bai_path = "/Users/yifeiwan/Projects/bamstorm/tests/chr1.bam.bai";
        let bai_path = "/Users/yifeiwan/Projects/bamstorm_old/test.bam.bai";
        let bam_path = "/Users/yifeiwan/Projects/bamstorm_old/test.bam";
        let bam_size = get_file_size(bam_path)?;
        let linear_indexes_all = get_linear_indexes(bai_path, bam_size)?;
        let intervals = get_linear_intervals(&linear_indexes_all)?;
        // let ind = bai::fs::read(bai_path)?;
        // let ref_seqs = ind.reference_sequences();
        // let ref_seq_0 = ref_seqs.iter().next();
        // println!("{:?}", ref_seq_0);
        assert_eq!(linear_indexes_all.len(), 176783);
        assert_eq!(intervals.len(), 176782);
        Ok(())
    }
}
