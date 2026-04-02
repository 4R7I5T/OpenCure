import json

# Define the file path for the uploaded jsonl file
file_path = "./OpenCureV1.2.jsonl"

# Define a function to clean sequences
def clean_sequence(sequence):
    # Remove invalid characters (like 'N') and ensure it's uppercase
    return ''.join([base for base in sequence if base in 'ACGT']).upper()

# Define a function to extract sequences from the jsonl file
def extract_sequences_from_jsonl(file_path):
    sequences = []
    with open(file_path, 'r') as jsonl_file:
        for line in jsonl_file:
            entry = json.loads(line.strip())
            # Assuming the file contains DNA sequences under specific keys, adjust as needed
            target_sequence = entry.get("target_sequence", "")
            guide_rna = entry.get("guide_rna", "")
            # Clean the sequences to remove invalid characters like 'N'
            clean_target_sequence = clean_sequence(target_sequence)
            clean_guide_rna = clean_sequence(guide_rna)
            if clean_target_sequence and clean_guide_rna:
                sequences.append((clean_target_sequence, clean_guide_rna))
    return sequences

# Define a function to create Cas-OFFinder input file
def create_cas_offinder_input(sequences, output_file):
    with open(output_file, 'w') as outfile:
        # Header for Cas-OFFinder format: Number of mismatches
        outfile.write("Nucleotide\tDNA sequence\tRNA sequence\n")
        for target_sequence, guide_rna in sequences:
            # Write each target sequence and guide RNA in the correct format
            outfile.write(f"DNA\t{target_sequence}\t{guide_rna}\n")

# Main function to extract sequences and create the Cas-OFFinder input file
def main():
    sequences = extract_sequences_from_jsonl(file_path)
    if sequences:
        output_file = "/mnt/data/cas_offinder_input.txt"
        create_cas_offinder_input(sequences, output_file)
        print(f"Cas-OFFinder input file created successfully at {output_file}")
    else:
        print("No valid sequences found in the JSONL file.")

# Run the main function
if __name__ == "__main__":
    main()