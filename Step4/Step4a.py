import requests
from Bio import Entrez, SeqIO
import pandas as pd
import io
import time
import subprocess
import os
import logging
from http.client import IncompleteRead

# Setup logging
logging.basicConfig(filename='process.log', level=logging.INFO, format='%(asctime)s %(message)s')

# Define your email and API key for NCBI
Entrez.email = "your_registered_email"
Entrez.api_key = "your_api_key"

# Function to fetch sequence from NCBI with retries
def fetch_sequence(chromosome, start, end, retries=3):
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(db="nucleotide", id=chromosome, seq_start=start, seq_end=end, rettype="fasta", retmode="text", api_key=Entrez.api_key)
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return str(record.seq)
        except IncompleteRead as e:
            logging.error(f"Incomplete read error fetching sequence for {chromosome}:{start}-{end}: {e}. Retrying...")
            attempt += 1
            time.sleep(5)  # Wait before retrying
        except Exception as e:
            logging.error(f"Error fetching sequence for {chromosome}:{start}-{end}: {e}. Retrying...")
            attempt += 1
            time.sleep(5)  # Wait before retrying
    return None

# Function to clean the fetched sequence
def clean_sequence(sequence):
    return ''.join(sequence.split()).upper()  # Remove whitespace and convert to uppercase

# Function to write sequences to a FASTA file for FlashFry
def write_fasta(sequences, filename="sequences.fasta"):
    with open(filename, 'w') as fasta_file:
        for i, seq in enumerate(sequences):
            if seq:  # Ensure sequence is not None
                fasta_file.write(f">sequence_{i}\n{seq}\n")

# Function to run FlashFry
def run_flashfry(input_fasta, output_dir="flashfry_output"):
    flashfry_path = "./FlashFry-assembly-1.15.jar"  # Update with the actual path to FlashFry
    genome_index = "./GRCh38_index.db"  # Update with the actual path to the genome index
    os.makedirs(output_dir, exist_ok=True)
    command = [
        "java", "-jar", flashfry_path,
        "discover",
        "--fasta", input_fasta,
        "--genome", genome_index,
        "--positionOutput", os.path.join(output_dir, "./results.tsv"),
        "--scoringMetrics", "Doench2016"
    ]
    subprocess.run(command, check=True)

# Function to parse FlashFry results
def parse_flashfry_results(results_file="./results.tsv"):
    gRNAs = {}
    if os.path.exists(results_file):
        with open(results_file, 'r') as f:
            for line in f:
                if line.startswith("sequence"):
                    continue
                parts = line.strip().split('\t')
                seq_id = parts[0]
                gRNA_seq = parts[1]
                gRNAs[seq_id] = gRNA_seq
    return gRNAs

# Load the CSV file with error handling for inconsistent fields
file_path = './resp003.csv'

# Read the file content directly using a context manager
with open(file_path, 'r') as file:
    file_content = file.readlines()

# Identify lines with inconsistent number of fields
field_counts = [len(line.split(',')) for line in file_content]
expected_fields = max(set(field_counts), key=field_counts.count)

# Filter out lines that do not match the expected number of fields
corrected_lines = [line for line in file_content if len(line.split(',')) == expected_fields]

# Create a corrected CSV data string
corrected_csv_data = '\n'.join(corrected_lines)

# Use StringIO to read the corrected CSV data into a pandas DataFrame
corrected_csv_data_io = io.StringIO(corrected_csv_data)
data = pd.read_csv(corrected_csv_data_io)

# Strip leading and trailing spaces from column names
data.columns = data.columns.str.strip()

# Verify the cleaned column names
logging.info("Column names: %s", data.columns)

# Initialize lists for sequences and results
sequences = []
sequence_ids = []
checkpoint_interval = 10
checkpoint_file = 'checkpoint.csv'

# Load checkpoint if it exists
if os.path.exists(checkpoint_file):
    checkpoint_data = pd.read_csv(checkpoint_file)
    sequences = checkpoint_data['Sequence'].tolist()
    sequence_ids = checkpoint_data['Sequence_ID'].tolist()
    processed_indices = set(checkpoint_data['Index'].tolist())
else:
    processed_indices = set()

# Function to save a checkpoint
def save_checkpoint():
    checkpoint_df = pd.DataFrame({
        'Index': list(processed_indices),
        'Sequence': sequences,
        'Sequence_ID': sequence_ids
    })
    checkpoint_df.to_csv(checkpoint_file, index=False)
    logging.info("Checkpoint saved.")

# Function to process each row sequentially
def process_row(index, row):
    if index in processed_indices:
        return None
    chromosome = row['Chromosome'].strip()
    position = int(row['Position'])
    variant_type = row['Variant_Type']
    reference = row['Reference']
    alternate = row['Alternate']
    
    # Define the region around the variant (125 bases upstream and downstream)
    start = max(1, position - 125)
    end = position + 125
    
    logging.info(f"Fetching sequence for {chromosome} from {start} to {end}")
    
    # Fetch the sequence
    sequence = fetch_sequence(chromosome, start, end)
    if not sequence:
        logging.error(f"Failed to fetch sequence for {chromosome} from {start} to {end}")
        return (index, None, f"sequence_{index}")
    
    # Clean the fetched sequence
    cleaned_sequence = clean_sequence(sequence)
    logging.info(f"Cleaned sequence for {chromosome} from {start} to {end}")
    
    return (index, cleaned_sequence, f"sequence_{index}")

# Process each row sequentially
for index, row in data.iterrows():
    result = process_row(index, row)
    if result:
        index, cleaned_sequence, seq_id = result
        sequences.append(cleaned_sequence)
        sequence_ids.append(seq_id)
        processed_indices.add(index)
        if len(processed_indices) % checkpoint_interval == 0:
            save_checkpoint()

# Write sequences to a FASTA file
write_fasta(sequences)
logging.info("Sequences written to FASTA file.")

# Run FlashFry
try:
    run_flashfry("sequences.fasta")
    logging.info("FlashFry executed successfully.")
except subprocess.CalledProcessError as e:
    logging.error(f"Error running FlashFry: {e}")

# Parse FlashFry results
gRNAs = parse_flashfry_results()

# Update the DataFrame with FlashFry results
crispor_input = []
crispor_output = []

for seq_id in sequence_ids:
    crispor_input.append(sequences[int(seq_id.split('_')[1])])
    crispor_output.append(gRNAs.get(seq_id, None))

# Add FlashFry input and output to the DataFrame
data['FlashFry_Input'] = crispor_input
data['FlashFry_Output'] = crispor_output

# Save the updated DataFrame to a new CSV file
updated_file_path = './updated_file.csv'
data.to_csv(updated_file_path, index=False)
logging.info("Updated DataFrame saved to CSV file.")

import ace_tools as tools; tools.display_dataframe_to_user(name="Updated Variant Data with gRNA Design", dataframe=data)