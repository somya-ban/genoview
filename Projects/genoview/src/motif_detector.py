# src/motif_detector.py

import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import traceback # For error handling in main block

# --- Define Real Motif Patterns (Using examples from before) ---
# Consider researching specific patterns relevant to your interests later.
MOTIF_PATTERNS = {
    # Example 1: Simplified TATA Box
    "TATA_BOX_LIKE": r"TATA[AT]A[AT]", # Using raw string prefix 'r'

    # Example 2: Basic GC Box
    "GC_BOX_LIKE": r"GGGCGG",

    # Example 3: Initiator Element (common alternative promoter element)
    # Consensus varies, e.g., YYANWYY (Y=[CT], N=[AGCT], W=[AT])
    "INR_LIKE": r"[CT][CT]A[AGCT][AT][CT][CT]",

    # Example 4: Basic CAAT Box (another promoter element)
    "CAAT_BOX_LIKE": r"GG[CT]CAATCT" # Example consensus
}

def find_motifs(sequence_record: SeqRecord, motif_patterns: dict = MOTIF_PATTERNS) -> list[dict]:
    """
    Finds predefined motifs in a DNA sequence using regular expressions.

    Searches both forward and reverse complement strands.

    Args:
        sequence_record (SeqRecord): Biopython SeqRecord containing the DNA sequence.
        motif_patterns (dict): Dictionary where keys are motif IDs/names and values
                               are Python regex strings.

    Returns:
        list[dict]: List of found motifs with 'motif_id', 'start', 'end',
                    'strand', 'matched_sequence'. Coordinates are 1-based.
                    Returns empty list if none found or on error.
    """
    results = []
    try:
        seq_fwd = str(sequence_record.seq).upper() # Work with uppercase for consistency
        seq_rev = str(sequence_record.seq.reverse_complement()).upper()
        seq_len = len(seq_fwd)

        if seq_len == 0:
            print("Warning: Input sequence is empty. Cannot find motifs.")
            return []

        for motif_id, pattern in motif_patterns.items():
            # Compile regex for efficiency and error checking
            try:
                # Case-insensitive search usually makes sense for DNA motifs
                compiled_pattern = re.compile(pattern, re.IGNORECASE)
            except re.error as e:
                print(f"Warning: Invalid regex pattern for motif '{motif_id}': {e}. Skipping.")
                continue # Skip this pattern

            # Search forward strand
            for match in compiled_pattern.finditer(seq_fwd):
                start_0based = match.start()
                end_0based = match.end() # Exclusive end in regex match
                results.append({
                    "motif_id": motif_id,
                    "start": start_0based + 1, # Convert to 1-based start
                    "end": end_0based,         # 1-based inclusive end
                    "strand": '+',
                    "matched_sequence": match.group(0)
                })

            # Search reverse strand
            for match in compiled_pattern.finditer(seq_rev):
                rev_start_0based = match.start()
                rev_end_0based = match.end() # Exclusive index

                # Map to 1-based inclusive coordinates on forward strand
                # Start on Fwd = Length - End_on_Rev + 1
                # End on Fwd   = Length - Start_on_Rev
                fwd_start_1based = seq_len - rev_end_0based + 1
                fwd_end_1based = seq_len - rev_start_0based

                results.append({
                    "motif_id": motif_id,
                    "start": fwd_start_1based,
                    "end": fwd_end_1based,
                    "strand": '-',
                    # Store the sequence matched on the REVERSE COMPLEMENT strand
                    "matched_sequence": match.group(0)
                    # Optional: Store fwd strand sequence? sequence_record.seq[fwd_start_1based-1:fwd_end_1based]
                })

        # Sort results by start position for readability
        results.sort(key=lambda x: x['start'])
        return results

    except Exception as e:
        print(f"An error occurred during motif finding:")
        traceback.print_exc()
        return [] # Return empty list on error

# --- Example Usage ---
if __name__ == "__main__":
    # Create a test sequence containing some motifs on both strands
    # Fwd: ...TATAATAA...GGGCGG... CCAATCT...
    # Rev: ...AGATTGG...CCGCC...TTATTATA...
    test_dna = ("AGCTTAGC" + "TATAATAA" + # TATA_BOX_LIKE (+) @ 9-16
                "GCTAGCTACG" + "GGGCGG" +   # GC_BOX_LIKE (+) @ 27-32
                "TAGCTAGCT" + "CCGCCCTA" + # Contains complement of GC_BOX_LIKE (-) -> CCGCCC @ 41 -> maps to 43-48
                "GCTAG" + "GGCCAATCT" + # CAAT_BOX_LIKE (+) @ 54-63
                "AGCT") # Total length 67
    # Let's check reverse mapping for GC_BOX_LIKE: CCGCCC found on FWD at 41-46.
    # RevComp: AGCT + AGATTGGCC + CTAGC + TAGGGCGG + CTAGCTAGC + TTATTATA + GCTAAGCT (Len 67)
    # Search for GGGCGG in RevComp: Found at index 21-26 (0-based)
    # fwd_start = 67 - 27 + 1 = 41
    # fwd_end   = 67 - 21 = 46. Correct!

    test_record = SeqRecord(Seq(test_dna), id="test_motif_seq")

    print(f"Finding motifs in: {test_record.id} (Length: {len(test_record.seq)})")
    # print(f"Sequence: {test_record.seq}")
    # print(f"RevComp : {test_record.seq.reverse_complement()}")

    found_motifs = find_motifs(test_record) # Use default patterns

    if found_motifs:
        print(f"\nFound {len(found_motifs)} motif instance(s):")
        for motif in found_motifs:
            print(f" - ID: {motif['motif_id']}, Start: {motif['start']}, End: {motif['end']}, Strand: {motif['strand']}, Match: {motif['matched_sequence']}")
            # Verify the matched sequence on the forward strand for '-' strand hits
            if motif['strand'] == '-':
                original_seq_slice = test_record.seq[motif['start']-1:motif['end']]
                print(f"     (Fwd strand: {original_seq_slice} -> RevComp: {original_seq_slice.reverse_complement()})") # Check if RevComp matches the 'matched_sequence'
    else:
        print("\nNo predefined motifs found.")