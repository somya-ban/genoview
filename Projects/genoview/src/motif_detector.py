# src/motif_detector.py
"""
Detects predefined sequence motifs using regular expressions.

This module provides a function to scan a DNA sequence (both strands)
for occurrences of known motifs defined by regular expression patterns.
"""

import re
import traceback
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from typing import List, Dict, Any # For type hinting

# --- Default Motif Patterns ---
# Define biologically relevant motifs using Python regex syntax.
# Keys are descriptive names/IDs, values are raw regex strings (r"...")
# Using IGNORECASE search, so patterns can be defined in uppercase.
# Examples:
# TATA Box: Common core promoter element
# GC Box: Promoter element, often binds Sp1 transcription factor
# CAAT Box: Promoter element, influences transcription rate
# INR Element: Initiator element, alternative to TATA box
MOTIF_PATTERNS = {
    "TATA_BOX_LIKE": r"TATA[AT]A[AT]",
    "GC_BOX_LIKE": r"GGGCGG",
    "INR_LIKE": r"[CT][CT]A[AGCT][AT][CT][CT]", # Example consensus YYANWYY
    "CAAT_BOX_LIKE": r"GG[CT]CAATCT" # Example consensus
    # Add more motifs here if desired
}

def find_motifs(sequence_record: SeqRecord, motif_patterns: Dict[str, str] = MOTIF_PATTERNS) -> List[Dict[str, Any]]:
    """
    Finds occurrences of predefined motifs in a DNA sequence using regex.

    Scans both the forward (+) and reverse complement (-) strands of the input
    sequence for matches against the provided regex patterns. Coordinates are
    reported relative to the forward strand (1-based).

    Args:
        sequence_record: Biopython SeqRecord containing the DNA sequence.
        motif_patterns: A dictionary where keys are motif IDs/names (str)
                        and values are Python regular expression strings (str)
                        representing the patterns to search for. Defaults to
                        the globally defined `MOTIF_PATTERNS`.

    Returns:
        A list of dictionaries, where each dictionary represents a found
        motif instance with the following keys:
        - 'motif_id' (str): The ID/name of the matched motif pattern.
        - 'start' (int): The 1-based start coordinate on the forward strand.
        - 'end' (int): The 1-based inclusive end coordinate on the forward strand.
        - 'strand' (str): '+' if found on the forward strand, '-' if found
                          on the reverse complement strand.
        - 'matched_sequence' (str): The actual sequence segment that matched
                                     the pattern (from the strand it was found on).

        Returns an empty list if no motifs are found or if an error occurs.
    """
    results: List[Dict[str, Any]] = []
    try:
        # Prepare sequences: work with uppercase for case-insensitive regex matching
        seq_fwd: str = str(sequence_record.seq).upper()
        seq_rev: str = str(sequence_record.seq.reverse_complement()).upper()
        seq_len: int = len(seq_fwd)

        if seq_len == 0:
            print("Warning: Input sequence is empty. Cannot find motifs.")
            return []

        # Iterate through each motif pattern provided
        for motif_id, pattern in motif_patterns.items():
            # Compile the regex pattern for efficiency and validation.
            # Use IGNORECASE for case-insensitive matching.
            try:
                compiled_pattern: re.Pattern = re.compile(pattern, re.IGNORECASE)
            except re.error as e:
                # Warn if a pattern is invalid and skip it
                print(f"Warning: Invalid regex pattern for motif '{motif_id}': {e}. Skipping.")
                continue

            # --- Search Forward Strand (+) ---
            # Use finditer to get match objects with start/end positions
            for match in compiled_pattern.finditer(seq_fwd):
                start_0based: int = match.start()
                end_0based: int = match.end() # Regex end() is exclusive index

                results.append({
                    "motif_id": motif_id,
                    "start": start_0based + 1, # Convert to 1-based start
                    "end": end_0based,         # End is 1-based inclusive
                    "strand": '+',
                    "matched_sequence": match.group(0) # The actual matched string
                })

            # --- Search Reverse Complement Strand (-) ---
            for match in compiled_pattern.finditer(seq_rev):
                rev_start_0based: int = match.start()
                rev_end_0based: int = match.end() # Exclusive index on reverse string

                # Map coordinates back to the original forward strand (1-based inclusive)
                # Start on Fwd = Seq Length - End_on_Rev (exclusive 0-based) + 1
                fwd_start_1based: int = seq_len - rev_end_0based + 1
                # End on Fwd   = Seq Length - Start_on_Rev (inclusive 0-based)
                fwd_end_1based: int = seq_len - rev_start_0based

                results.append({
                    "motif_id": motif_id,
                    "start": fwd_start_1based,
                    "end": fwd_end_1based,
                    "strand": '-',
                    # Store the sequence segment matched on the REVERSE COMPLEMENT strand
                    "matched_sequence": match.group(0)
                })

        # Sort final results by start position for better readability
        results.sort(key=lambda x: x['start'])
        return results

    except Exception as e:
        # Catch-all for any unexpected errors during the process
        print(f"An error occurred during motif finding:")
        traceback.print_exc()
        return [] # Return empty list on error

# --- Example Usage Block ---
if __name__ == "__main__":
    print("--- Testing motif_detector.py ---")

    # Test sequence designed to contain examples of default motifs
    test_dna = ("AGCTTAGC" + "TATAATAA" + # TATA_BOX_LIKE (+) @ 9-16
                "GCTAGCTACG" + "GGGCGG" +   # GC_BOX_LIKE (+) @ 27-32
                "TAGCTAGCT" + "CCGCCCTA" + # Contains complement of GC_BOX_LIKE (-) -> GGGCGG @ 41-46
                "GCTAG" + "GGCCAATCT" + # CAAT_BOX_LIKE (+) @ 54-63
                "AGCT") # Total length 67
    test_record = SeqRecord(Seq(test_dna), id="test_motif_seq")

    print(f"\nTesting sequence: {test_record.id} (Length: {len(test_record.seq)})")
    print(f"Sequence: {test_record.seq}")
    print(f"RevComp : {test_record.seq.reverse_complement()}")

    # Call the function with default patterns
    found_motifs = find_motifs(test_record)

    if found_motifs:
        print(f"\nFound {len(found_motifs)} motif instance(s):")
        for motif in found_motifs:
            print(f" - ID: {motif['motif_id']}, Start: {motif['start']}, End: {motif['end']}, "
                  f"Strand: {motif['strand']}, Match: {motif['matched_sequence']}")
            # For reverse strand hits, verify the mapping by showing the forward strand segment
            if motif['strand'] == '-':
                try:
                    original_seq_slice = test_record.seq[motif['start']-1:motif['end']]
                    print(f"     (Fwd strand: {original_seq_slice} -> RevComp check: {original_seq_slice.reverse_complement()})")
                except IndexError:
                     print("     (Error slicing original sequence for verification)")
    else:
        print("\nNo predefined motifs found in the test sequence.")