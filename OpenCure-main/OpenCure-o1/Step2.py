import tiktoken
import os

def count_tokens_in_file(file_path, encoding_name="cl100k_base", chunk_size=1024 * 1024):
    """
    Counts the number of tokens in a given file.

    Args:
        file_path (str): Path to the file.
        encoding_name (str): Name of the encoding to use for tokenization. Default is "cl100k_base" for GPT-4.
        chunk_size (int): Size of each chunk to read from the file. Default is 1MB.

    Returns:
        int: Total number of tokens in the file.
    """
    # Check if the file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    # Initialize the tokenizer
    tokenizer = tiktoken.get_encoding(encoding_name)

    total_tokens = 0

    # Read the file in chunks to handle large files
    with open(file_path, 'r', encoding='utf-8') as file:
        while True:
            chunk = file.read(chunk_size)
            if not chunk:
                break
            tokens = tokenizer.encode(chunk)
            total_tokens += len(tokens)

    return total_tokens

if __name__ == "__main__":
    file_path = 'prompt01.txt'  # Use the file generated in step1.py
    try:
        total_tokens = count_tokens_in_file(file_path)
        print(f"Number of tokens in the file '{file_path}': {total_tokens}")
    except FileNotFoundError as e:
        print(e)
    except Exception as e:
        print(f"An error occurred: {e}")