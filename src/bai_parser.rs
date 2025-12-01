use noodles::bam::bai;
use noodles::bgzf::VirtualPosition;

// Get all linear indexes from a BAI file
pub fn get_linear_indexes(bai_path: &str) -> Vec<VirtualPosition> {
    // Read the BAI file
    let bai_file = bai::fs::read(bai_path).expect("Failed to read BAI file");

    // Prepare a vector to store all linear indexes
    let mut linear_indexes = Vec::new();

    // Iterate over all reference sequences in the BAI
    for ref_seq in bai_file.reference_sequences() {
        // Clone the linear index for this reference sequence
        let mut linear_index = ref_seq.index().clone();

        // Remove duplicate entries
        linear_index.dedup();

        // Append to the master list
        linear_indexes.extend(linear_index);
    }

    linear_indexes.sort();
    linear_indexes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_bai() -> Result<(), Box<dyn std::error::Error>> {
        let bai_path = "E:/project/bamstrom/tests/mt.sorted.bam.bai";
        let linear_indexes_all = get_linear_indexes(bai_path);
        // let ind = bai::fs::read(bai_path)?;
        // let ref_seqs = ind.reference_sequences();
        // let ref_seq_0 = ref_seqs.iter().next();
        // println!("{:?}", ref_seq_0);
        println!("All indexes >>> {:?}", linear_indexes_all.len());
        Ok(())
    }
}