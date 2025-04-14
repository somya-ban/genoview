# src/fasta_parser.py
"""
Parses FASTA formatted sequence files.

This module provides functionality to read and validate sequences from files
adhering to the FASTA format specification. It uses the Biopython library
for robust parsing.
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import Iterator, Optional # For type hinting

def parse_fasta_file(filepath: str) -> Iterator[SeqRecord]:
    """
    Parses a FASTA file and yields SeqRecord objects for each sequence found.

    This function handles file existence checks and basic FASTA format validation
    by attempting to parse sequences using Biopython's SeqIO module. It yields
    records one by one, making it memory-efficient for large files.

    Args:
        filepath: The absolute or relative path to the input FASTA file.

    Yields:
        Bio.SeqRecord.SeqRecord: A SeqRecord object for each sequence entry
                                  discovered in the FASTA file.

    Raises:
        FileNotFoundError: If the specified file does not exist at `filepath`.
        ValueError: If the file exists but is empty, contains no valid FASTA
                    sequences, or if Biopython encounters a parsing error
                    suggesting an invalid format.
        IOError: If there's a general issue reading the file.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Error: Input FASTA file not found at '{filepath}'")

    # Use a try-except block to catch potential issues during file handling/parsing
    try:
        # SeqIO.parse returns an iterator, efficient for large files.
        records_iterator = SeqIO.parse(filepath, "fasta")

        # Check if the iterator yields at least one record to ensure the file
        # isn't empty or completely malformed.
        first_record: Optional[SeqRecord] = None
        try:
            first_record = next(records_iterator)
        except StopIteration:
            # This means the file was opened but contained no sequences (or only comments)
            raise ValueError(f"Error: No sequences found in '{filepath}'. File might be empty or invalid.")

        # If we successfully got the first record, yield it first.
        yield first_record

        # Then, yield the rest of the records from the iterator.
        yield from records_iterator

    except ValueError as e:
        # Re-raise ValueErrors (like the one from StopIteration above or from Biopython)
        # with more context. Biopython might raise ValueError for format issues.
        raise ValueError(f"Error parsing FASTA file '{filepath}': {e}. Please ensure it's a valid FASTA format.") from e
    except IOError as e:
        # Catch file reading errors
        raise IOError(f"Error reading file '{filepath}': {e}") from e
    except Exception as e:
        # Catch any other unexpected errors during parsing
        raise RuntimeError(f"An unexpected error occurred while parsing '{filepath}': {e}") from e

# --- Example Usage Block ---
# This block demonstrates how to use the parser function when the script
# is executed directly. It's useful for basic testing of the parser itself.
if __name__ == "__main__":
    import sys # Needed for traceback printing

    print("--- Testing fasta_parser.py ---")

    # Assumes running from the project root directory (e.g., 'genoview/')
    # Use a sample file expected to be in the 'data/' subdirectory
    test_fasta_path = os.path.join("data", "sample.fasta")
    print(f"\nAttempting to parse: {test_fasta_path}")

    try:
        sequence_count = 0
        # Iterate through the SeqRecord objects yielded by the parser
        for record in parse_fasta_file(test_fasta_path):
            sequence_count += 1
            print("-" * 20)
            print(f"Sequence {sequence_count}:")
            print(f"  ID: {record.id}")
            print(f"  Description: {record.description}")
            print(f"  Length: {len(record.seq)} bp")
            # Display only the start of the sequence to avoid large outputs
            print(f"  Sequence (first 60 bp): {record.seq[:60]}...")

        print("-" * 20)
        if sequence_count > 0:
            print(f"\nSuccessfully parsed {sequence_count} sequence(s).")
        else:
            # This path should ideally not be reached if the function raises errors correctly
            print("\nParsing finished, but no sequences were yielded (check error messages above).")

    except (FileNotFoundError, ValueError, IOError, RuntimeError) as e:
        print(f"\nERROR during parsing test: {e}")
    except Exception as e:
        # Catch any truly unexpected errors during the test run itself
        print(f"\nUNEXPECTED ERROR during testing:")
        traceback.print_exc() # Print full traceback for debugging unexpected issues

    # Example of testing a non-existent file
    print("\nAttempting to parse non-existent file...")
    try:
        for record in parse_fasta_file("non_existent_file.fasta"):
            pass # Should not reach here
    except FileNotFoundError as e:
        print(f"Successfully caught expected error: {e}")
    except Exception as e:
        print(f"Caught unexpected error instead of FileNotFoundError: {e}")