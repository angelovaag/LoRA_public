import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Sum Reads based on TAXid for each sample")
    parser.add_argument("-d", "--input_dir", required=True, help="Path to the directory containing input files (_4krona.txt [contigName Reads TAXid])")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output file (merged table)")
    return parser.parse_args()

def main():
    args = parse_args()

    # List all files in the specified directory
    input_files = [f for f in os.listdir(args.input_dir) if f.endswith(".txt")]

    # Dictionary to store the sum of Reads for each TAXid and sample
    taxid_reads = {}

    # Process each input file
    for input_file in input_files:
         # Extract sample name from the file name and remove "_4krona"
        sample_name = os.path.splitext(input_file)[0].replace("_4krona", "")
        
        with open(os.path.join(args.input_dir, input_file)) as file:
            next(file)  # Skip the header line
            for line in file:
                contig_name, reads, taxid = line.strip().split("\t")
                key = (taxid, sample_name)
                taxid_reads[key] = taxid_reads.get(key, 0) + int(reads)

    # Get unique TAXids and sample names
    unique_taxids = sorted(set(taxid for taxid, _ in taxid_reads.keys()))
    unique_samples = sorted(set(sample for _, sample in taxid_reads.keys()))

    # Write the output to the specified file
    with open(args.output_file, "w") as output_file:
        # Print the header
        output_file.write("\t".join(["TAXid"] + unique_samples) + "\n")

        # Print the data
        for taxid in unique_taxids:
            values = [str(taxid_reads.get((taxid, sample), 0)) for sample in unique_samples]
            output_file.write("\t".join([taxid] + values) + "\n")

if __name__ == "__main__":
    main()
