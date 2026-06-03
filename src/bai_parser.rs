use std::fs::File;

use noodles::bam::bai;
use noodles::bgzf::VirtualPosition;

const DEFALUT_MAX_COMPRESSED_OFFSET: u64 = (1 << 48) - 1;
// const DEFALUT_MIN_UNCOMPRESSED_OFFSET: u16 = u16::MAX;

// Get file size
pub fn get_file_size(file_path: &str) -> std::io::Result<u64> {
    let file = File::open(file_path)?;
    Ok(file.metadata().unwrap().len())
}

// Get all linear indexes from a BAI file
pub fn get_linear_indexes(
    bai_path: &str,
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
    if linear_indexes.last().unwrap().compressed() > DEFALUT_MAX_COMPRESSED_OFFSET {
        return Err(format!(
            "The max coffset {:?} exceeds MAX coffset {:?}",
            linear_indexes.last().unwrap().compressed(),
            DEFALUT_MAX_COMPRESSED_OFFSET
        )
        .into());
    }
    // Create a VirtualPosition at the very start of the file
    // let start_voffset = VirtualPosition::new(0, 0).unwrap();

    // Insert at the beginning
    // linear_indexes.insert(0, start_voffset);

    Ok(linear_indexes)
}

pub fn reduce_linear_indexes(linear_indexes: &[VirtualPosition]) -> Vec<VirtualPosition> {
    let mut reduced: Vec<VirtualPosition> = Vec::with_capacity(linear_indexes.len());

    for &voffset in linear_indexes {
        match reduced.last_mut() {
            Some(last) if last.compressed() == voffset.compressed() => {
                // Same compressed offset, keep the smaller virtual offset
                if voffset < *last {
                    *last = voffset;
                }
            }
            _ => {
                reduced.push(voffset);
            }
        }
    }

    reduced
}

// TODO: Some intervals have same start coffsets, need to combine them or reuse the block somehow
pub fn get_linear_intervals(
    linear_indexes: &[VirtualPosition],
) -> Result<Vec<(VirtualPosition, VirtualPosition)>, String> {
    // Ensure there are enough indexes to form intervals
    if linear_indexes.len() < 2 {
        return Err("Not enough linear indexes to form intervals".to_string());
    }
    // Create intervals from consecutive linear indexes
    let reduced_indexes = reduce_linear_indexes(&linear_indexes);
    let mut intervals = Vec::new();
    for window in reduced_indexes.windows(2) {
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
    #[ignore]
    fn test_read_bai() -> Result<(), Box<dyn std::error::Error>> {
        let bai_path = "/Users/yifeiwan/Projects/bamstorm/tests/chr1.bam.bai";
        // let bai_path = "/Users/yifeiwan/Projects/bamstorm_old/test.bam.bai";
        let linear_indexes_all = get_linear_indexes(bai_path)?;
        println!(
            "First linear index: {:?}",
            linear_indexes_all.first().unwrap()
        );
        let intervals = get_linear_intervals(&linear_indexes_all)?;
        println!("First interval: {:?}", intervals[0]);
        Ok(())
    }

    #[test]
    #[ignore]
    fn test_reduce_linear_indexes_merges_same_compressed_offsets() {
        // Input virtual offsets are already sorted
        let linear_indexes = vec![
            VirtualPosition::new(1, 0).unwrap(),
            VirtualPosition::new(1, 10).unwrap(),
            VirtualPosition::new(1, 20).unwrap(),
            VirtualPosition::new(3, 0).unwrap(),
            VirtualPosition::new(5, 5).unwrap(),
        ];

        let reduced = reduce_linear_indexes(&linear_indexes);

        // Expect only one entry per compressed offset
        assert_eq!(reduced.len(), 3);

        assert_eq!(reduced[0].compressed(), 1);
        assert_eq!(reduced[0].uncompressed(), 0);

        assert_eq!(reduced[1].compressed(), 3);
        assert_eq!(reduced[1].uncompressed(), 0);

        assert_eq!(reduced[2].compressed(), 5);
        assert_eq!(reduced[2].uncompressed(), 5);
    }

    #[test]
    #[ignore]
    fn check_get_linear_indexes() -> Result<(), Box<dyn std::error::Error>> {
        let bai_path = "/Users/yifeiwan/Projects/bamstorm/tests/chr1.bam.bai";
        let linear_indexes_all = get_linear_indexes(bai_path)?;
        let intervals = get_linear_intervals(&linear_indexes_all)?;
        println!("Head of linear indexes: {:?}", &intervals[0..2]);
        Ok(())
    }
}
