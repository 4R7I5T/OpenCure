import pandas as pd
import os

# Paths to your CSV files
input_file_path = "./output_30.csv"
output_dir = "./samples"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Specify the chunk size (number of lines to read at a time)
chunk_size = 10000

# Initialize a dictionary to hold the sequences data
sequence_data = {}

# Process the CSV file in chunks
for chunk in pd.read_csv(input_file_path, chunksize=chunk_size):
    for index, row in chunk.iterrows():
        sequence_id = row['chromosome']  # Replace with the correct column name for sequence ID
        on_site_score = row['on_site_score']  # Replace with the correct column name for on_site_score
        
        if sequence_id not in sequence_data:
            sequence_data[sequence_id] = []

        # Append the row as a tuple (on_site_score, row_data)
        sequence_data[sequence_id].append((on_site_score, row))
        
        # If we've collected more than 1000 rows for this sequence, sort by on_site_score and write the top 1000
        if len(sequence_data[sequence_id]) > 1000:
            sequence_data[sequence_id].sort(reverse=True, key=lambda x: x[0])
            top_data = sequence_data[sequence_id][:1000]
            sequence_data[sequence_id] = sequence_data[sequence_id][1000:]  # Keep the remaining for future processing

            # Write out the top data incrementally
            output_file_path = os.path.join(output_dir, f"{sequence_id}_top1000.csv")
            pd.DataFrame([row for _, row in top_data]).to_csv(output_file_path, mode='a', header=not os.path.exists(output_file_path), index=False)

# Write out any remaining data that didn't reach the 1000-row threshold
for sequence_id, data in sequence_data.items():
    if data:
        data.sort(reverse=True, key=lambda x: x[0])
        output_file_path = os.path.join(output_dir, f"{sequence_id}_top1000.csv")
        pd.DataFrame([row for _, row in data]).to_csv(output_file_path, mode='a', header=not os.path.exists(output_file_path), index=False)

print("Incremental sampling and saving completed.")