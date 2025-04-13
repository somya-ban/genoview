import os
from Bio import SeqIO # Import the SeqIO module

def parse_fasta_file(filepath):
    """
    Parses a FASTA file and yields SeqRecord objects.

    Args:
        filepath (str): The path to the FASTA file.

    Yields:
        Bio.SeqRecord.SeqRecord: A SeqRecord object for each sequence in the file.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file is not in a valid FASTA format.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Error: File not found at '{filepath}'")

    try:
        # Use SeqIO.parse for files with potentially multiple sequences
        # It returns an iterator, which is memory-efficient for large files
        records = SeqIO.parse(filepath, "fasta")
        # We need to check if the file actually contains any records
        # Converting iterator to list will load all records, use carefully for large files
        # For now, let's just get the first record to check format validity
        first_record = next(records, None) # Get first item or None if empty

        if first_record is None:
             # Handle empty file or file with only comments
             raise ValueError(f"Error: No sequences found in '{filepath}'. Is it a valid FASTA file?")

        # If the first record was valid, we can now iterate properly.
        # We need to chain the first record back with the rest of the iterator.
        # (Alternatively, re-open the file, but this is more efficient)
        from itertools import chain
        all_records = chain([first_record], records)

        # Yield each record
        yield from all_records

    except Exception as e:
        # Catch potential parsing errors from Biopython or other issues
        raise ValueError(f"Error parsing FASTA file '{filepath}': {e}. Please ensure it's a valid FASTA format.")

# --- Example Usage ---
if __name__ == "__main__":
    # This block runs only when the script is executed directly
    # Adjust the path if your script is run from a different directory
    # Assumes running from the main 'genoview' directory: python src/fasta_parser.py
    fasta_file_path = os.path.join("data", "sample.fasta") # Construct path relative to project root

    print(f"Attempting to parse FASTA file: {fasta_file_path}")

    try:
        sequence_count = 0
        for record in parse_fasta_file(fasta_file_path):
            sequence_count += 1
            print("-" * 20)
            print(f"Sequence ID: {record.id}")
            print(f"Description: {record.description}")
            print(f"Sequence Length: {len(record.seq)} bases")
            # Don't print the whole sequence for potentially long files
            print(f"Sequence (first 30 bases): {record.seq[:30]}...")
            # You can access the full sequence via record.seq

        print("-" * 20)
        print(f"\nSuccessfully parsed {sequence_count} sequence(s) from the file.")

    except FileNotFoundError as e:
        print(e)
    except ValueError as e:
        print(e)
    except Exception as e:
        # Catch any other unexpected errors
        print(f"An unexpected error occurred: {e}")