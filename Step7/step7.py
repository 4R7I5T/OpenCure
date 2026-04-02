import pandas as pd
import re
import json

def load_crispr_data(csv_file):
    try:
        # Load the CSV file
        crispr_df = pd.read_csv(csv_file)
        return crispr_df
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return None

def parse_fasta(fasta_file):
    try:
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
        sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
        return sequence
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return None

def generate_input_output_scenarios(crispr_df, sequence):
    scenarios = []
    
    for _, row in crispr_df.iterrows():
        try:
            start_pos = int(row['start_pos'])
            
            # Extract a 100-character snippet from the FASTA sequence based on the provided starting position
            fasta_snippet = sequence[start_pos-1:start_pos+99]
            
            # Create a scenario with the user providing a sequence snippet and the assistant providing the required fields
            scenario = {
                "messages": [
                    {"role": "user", "content": f"Sequence starting at position {start_pos}: {fasta_snippet}"},
                    {"role": "assistant", "content": f"CRISPR ID: {row['crispr_id']}, Sequence: {row['sequence']}, Long Sequence: {row['long_sequence']}, Start Position: {row['start_pos']}, End Position: {row['end_pos']}, Cutsite: {row['cutsite']}"}
                ]
            }
            scenarios.append(scenario)
        
        except ValueError as ve:
            print(f"Skipping entry due to conversion error: {ve}")
        except KeyError as ke:
            print(f"Skipping entry due to missing data: {ke}")
    
    return scenarios

def save_scenarios_to_jsonl(scenarios, output_file):
    with open(output_file, 'w', encoding='utf-8') as f:
        for scenario in scenarios:
            json_line = json.dumps(scenario)
            f.write(json_line + '\n')

def main():
    fasta_file = 'sequence_0_region.fa'  # Ensure the correct path to your FASTA file
    csv_file = 'Seq0_CRISPR_1.csv'  # Path to the cleaned CSV
    output_file = 'output.jsonl'
    
    # Load the CRISPR data
    crispr_df = load_crispr_data(csv_file)
    if crispr_df is None:
        print("No CRISPR data loaded.")
        return
    
    # Parse the sequence from the FASTA file
    sequence = parse_fasta(fasta_file)
    if not sequence:
        print("No sequence data loaded.")
        return
    
    # Generate input/output scenarios
    scenarios = generate_input_output_scenarios(crispr_df, sequence)
    
    # Save the scenarios to a .jsonl file
    if scenarios:
        save_scenarios_to_jsonl(scenarios, output_file)
        print(f"Output saved to {output_file}")
    else:
        print("No scenarios generated.")

if __name__ == '__main__':
    main()