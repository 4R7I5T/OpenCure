import os
from openai import OpenAI
import pandas as pd
import tiktoken
import re
import csv

# Initialize the client using the OPENAI_API_KEY environment variable
api_key = os.getenv("OPENAI_API_KEY")
if not api_key:
    raise RuntimeError("OPENAI_API_KEY environment variable is not set.")
client = OpenAI(api_key=api_key)

job = '000'
prompt_path = "./prompt01.txt"
prompt_outputfile = f"resp{job}.csv"

# Set up initial prompt
prompt = """
Please return the requested output per the question that follows in a zipped file. Here is the question:

"Based on the following annotated variants from a tumor-normal comparison, suggest specific CRISPR-Cas9 gene edits that can revert the tumor genome to the normal state":
"""

endprompt = """This completes all data chunks, please process and return the requested output to the question I asked. 
Your response should ONLY include the CSV formatted data with these columns 'Chromosome, Position, Variant_Type, 
Reference, Alternate, Suggested_Edit'.
"""

# Function to tokenize text and split into batches
def batch_data_by_tokens(data, tokenizer, max_tokens):
    tokens = tokenizer.encode(data)
    for i in range(0, len(tokens), max_tokens):
        yield tokenizer.decode(tokens[i:i + max_tokens])

# Read the input data
with open(prompt_path, 'r') as file:
    crispr_data = file.read()

# Initialize the tokenizer
tokenizer = tiktoken.encoding_for_model("gpt-4")

# Define the maximum number of tokens per batch
max_tokens_per_batch = 128000 - len(tokenizer.encode(prompt)) - len(tokenizer.encode(endprompt)) - 100  # buffer for safety

results = []

# Process each batch
for batch in batch_data_by_tokens(crispr_data, tokenizer, max_tokens_per_batch):
    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": f"{prompt}\n{batch}\n{endprompt}"}
        ]
    )

    # Extract the content from the response
    content = resp.choices[0].message.content
    
    # Remove the markdown markers ⁠ csv and  ⁠
    content = content.replace('⁠  csv\n', '').replace('  ⁠', '')

    # Split the content into lines and filter out any empty lines
    lines = [line for line in content.split('\n') if line.strip()]
    
    # Add the lines to the results
    results.extend(lines)

# Write the results to a CSV file with a single clean header
header = ["Chromosome", "Position", "Variant_Type", "Reference", "Alternate", "Suggested_Edit"]

with open(prompt_outputfile, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header)

    for line in results:
        line = line.strip()
        if not line:
            continue
        # Skip any header-like lines returned by the model
        if line.lower().startswith("chromosome"):
            continue
        parts = [field.strip() for field in line.split(',')]
        if len(parts) != 6:
            continue
        try:
            int(parts[1])
        except ValueError:
            continue
        csvwriter.writerow(parts)

print(f"Data has been written to {prompt_outputfile}")