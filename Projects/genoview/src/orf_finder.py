# src/orf_finder.py
"""
Finds Open Reading Frames (ORFs) within biological sequences.

This module implements ORF detection using standard genetic codes and start/stop
codons based purely on the Biopython library. It scans all six reading frames.
"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data import CodonTable
import traceback # For printing detailed errors in __main__
from typing import List, Dict, Any, Optional # For type hinting

def find_orfs_biopython(sequence_record: SeqRecord, min_protein_length: int = 25, table: int = 1) -> List[Dict[str, Any]]:
    """
    Finds Open Reading Frames (ORFs) in a DNA sequence using Biopython.

    Identifies potential protein-coding regions based on standard start ('M')
    and stop ('*') codons in all 6 reading frames (3 forward, 3 reverse).
    Only ORFs encoding proteins meeting the minimum length are returned.

    Args:
        sequence_record: A Biopython SeqRecord object containing the DNA sequence.
                         The sequence alphabet should ideally be DNA.
        min_protein_length: The minimum length (in amino acids) of the translated
                            protein (excluding the stop codon) required for an
                            ORF to be reported. Defaults to 25 aa.
        table: The NCBI Genetic Code table number (integer) to use for translation.
               Defaults to 1 (Standard Code). See Biopython documentation
               for available table numbers.

    Returns:
        A list of dictionaries. Each dictionary represents a found ORF and
        contains the following keys:
        - 'orf_id' (str): A generated unique ID for this ORF within the sequence
                          (e.g., "orf_1", "orf_2"). Added later by pipeline usually.
        - 'start' (int): The 1-based start coordinate of the ORF on the
                         original **forward** strand.
        - 'end' (int): The 1-based end coordinate of the ORF (inclusive) on the
                       original **forward** strand. Includes the stop codon position.
        - 'strand' (str): '+' for forward strand, '-' for reverse strand.
        - 'length_bp' (int): The length of the ORF in base pairs (including stop codon).
        - 'protein_sequence' (str): The translated amino acid sequence, excluding
                                    the stop codon.

        Returns an empty list if no ORFs meeting the criteria are found or if
        an error occurs during processing.
    """
    orf_results: List[Dict[str, Any]] = []
    seq: Seq = sequence_record.seq
    seq_len: int = len(seq)

    if seq_len < min_protein_length * 3 + 3: # Basic check: sequence too short?
        # Needs at least min_protein_length codons + 1 stop codon = (min+1)*3 bases
        # print("Sequence too short to contain required ORF length.")
        return []

    try:
        # Get the specified codon table from Biopython
        codon_table: CodonTable.CodonTable = CodonTable.ambiguous_dna_by_id[table]
    except KeyError:
        print(f"Warning: Codon table ID {table} not found. Using standard table 1.")
        try:
            codon_table = CodonTable.ambiguous_dna_by_id[1]
        except KeyError:
            # This should essentially never happen with table 1
            print("Error: Standard codon table 1 not found in Biopython. Cannot proceed.")
            return []

    # Define valid start codons based on the chosen table (DNA sequence)
    # Note: Biopython's translate method implicitly handles start codons,
    # we primarily check for 'M' in the translated protein sequence.
    # start_codons = codon_table.start_codons
    # stop_codons = codon_table.stop_codons

    # Iterate through both the forward sequence (+1) and its reverse complement (-1)
    for strand_value, nuc_seq in [(+1, seq), (-1, seq.reverse_complement())]:
        strand_char: str = '+' if strand_value == 1 else '-'

        # Iterate through the three possible reading frames (0, 1, 2)
        for frame in range(3):
            # Slice the sequence to start at the correct frame
            frame_seq: Seq = nuc_seq[frame:]
            frame_len: int = len(frame_seq)

            # Ensure there's enough sequence for at least one codon in this frame
            if frame_len < 3:
                continue

            protein_seq: str = ""
            # Store the 0-based start index (relative to nuc_seq) for each potential codon
            # This allows mapping protein positions back to DNA positions precisely.
            dna_indices: List[int] = []

            # Translate codon by codon to handle sequences not perfectly divisible by 3
            for i in range(0, frame_len - (frame_len % 3), 3):
                codon: Seq = frame_seq[i:i+3]
                # Record the starting DNA index for this codon within the current strand/frame view
                # The index is relative to the start of `nuc_seq` (which is either fwd or revcomp)
                dna_indices.append(frame + i)
                try:
                    # Translate using the specified table, standard stop symbol '*'
                    amino_acid: str = codon.translate(table=table, stop_symbol="*")
                    protein_seq += amino_acid
                except CodonTable.TranslationError:
                    # If a codon is ambiguous or invalid for the table, represent as '?'
                    # Alternatively, could use 'X' or skip. '?' highlights potential issues.
                    protein_seq += "?"

            # Skip if no protein sequence was generated (e.g., frame too short)
            if not protein_seq:
                continue

            # Scan the translated protein sequence for valid ORFs (M...*)
            prot_len: int = len(protein_seq)
            # Track the 0-based index in protein_seq where the current potential ORF started ('M')
            current_orf_start_prot_idx: int = -1

            for aa_index in range(prot_len):
                amino_acid_char: str = protein_seq[aa_index]

                # Found a start codon ('M') and not currently inside an ORF:
                if amino_acid_char == 'M' and current_orf_start_prot_idx == -1:
                    current_orf_start_prot_idx = aa_index

                # Found a stop codon ('*') while inside a potential ORF:
                elif amino_acid_char == '*' and current_orf_start_prot_idx != -1:
                    current_orf_end_prot_idx: int = aa_index
                    # Calculate protein length (excluding stop codon)
                    current_prot_orf_len: int = current_orf_end_prot_idx - current_orf_start_prot_idx

                    # Check if the protein meets the minimum length requirement
                    if current_prot_orf_len >= min_protein_length:
                        # Map protein indices back to DNA coordinates on the current strand (nuc_seq)

                        # Start DNA coord (1-based): (0-based DNA index of M codon) + 1
                        dna_orf_start_on_strand: int = dna_indices[current_orf_start_prot_idx] + 1
                        # End DNA coord (1-based, inclusive): (0-based DNA index of stop codon) + 3
                        dna_orf_end_on_strand: int = dna_indices[current_orf_end_prot_idx] + 3

                        # Extract the protein sequence for this ORF (M to *)
                        orf_protein: str = protein_seq[current_orf_start_prot_idx:current_orf_end_prot_idx]
                        # Calculate DNA length
                        orf_dna_len: int = dna_orf_end_on_strand - dna_orf_start_on_strand + 1

                        # Map coordinates back to the original forward strand if ORF was found on reverse complement
                        if strand_value == -1:
                            # Remember seq_len is the length of the original forward strand
                            final_dna_start: int = seq_len - dna_orf_end_on_strand + 1
                            final_dna_end: int = seq_len - dna_orf_start_on_strand + 1
                        else:
                            # Coordinates are already relative to the forward strand
                            final_dna_start: int = dna_orf_start_on_strand
                            final_dna_end: int = dna_orf_end_on_strand

                        # Store the results for this valid ORF
                        orf_results.append({
                            # 'orf_id' is added later by the pipeline for uniqueness
                            "start": final_dna_start,
                            "end": final_dna_end,
                            "strand": strand_char,
                            "length_bp": orf_dna_len,
                            "protein_sequence": orf_protein
                        })

                    # Reset ORF tracking after encountering a stop codon, regardless of length check
                    current_orf_start_prot_idx = -1

    # Sort the found ORFs by their start position on the forward strand
    orf_results.sort(key=lambda x: x['start'])
    return orf_results

# --- Example Usage Block ---
if __name__ == "__main__":
    print("--- Testing orf_finder.py ---")

    # --- Test Case 1: Forward ORF ---
    # Expected ORF: Frame 0, MKKKKKGPF (9 aa), DNA 13-42 (+) Len 30 bp
    test_dna_fwd = "GGGAAACCCTTTATGAAAAAAAAAAAAAAAGGGCCCTTTTAGGGCCC" # Len 47
    record_fwd = SeqRecord(Seq(test_dna_fwd), id="test_fwd_orf", description="Test for forward ORF")
    print(f"\nTesting sequence: {record_fwd.id} (Length: {len(record_fwd.seq)} bp)")
    min_len_fwd = 9 # Set minimum length to match expected ORF

    try:
        found_orfs_fwd = find_orfs_biopython(record_fwd, min_protein_length=min_len_fwd)
        if found_orfs_fwd:
            print(f"Found {len(found_orfs_fwd)} ORF(s) [min_len={min_len_fwd} aa]:")
            for i, orf in enumerate(found_orfs_fwd):
                print(f"  ORF {i+1}: Start={orf['start']}, End={orf['end']}, Strand={orf['strand']}, "
                      f"Len={orf['length_bp']}bp, ProtLen={len(orf['protein_sequence'])}aa, Prot={orf['protein_sequence']}")
        else:
            print(f"No ORFs found [min_len={min_len_fwd} aa].")
    except Exception as e:
        print(f"\nERROR during forward ORF test:")
        traceback.print_exc()

    # --- Test Case 2: Reverse ORF (or lack thereof) ---
    # Expected: No significant ORFs meeting length criteria based on previous debug
    test_dna_rev = "CCCGGGCTACATTTTTTTTTTTCATAAAGGCCCGGG" # Len 36
    record_rev = SeqRecord(Seq(test_dna_rev), id="test_rev_orf", description="Test for reverse ORF")
    print(f"\nTesting sequence: {record_rev.id} (Length: {len(record_rev.seq)} bp)")
    min_len_rev = 5 # Set minimum length for this test

    try:
        found_orfs_rev = find_orfs_biopython(record_rev, min_protein_length=min_len_rev)
        if found_orfs_rev:
             print(f"Found {len(found_orfs_rev)} ORF(s) [min_len={min_len_rev} aa]:")
             for i, orf in enumerate(found_orfs_rev):
                 print(f"  ORF {i+1}: Start={orf['start']}, End={orf['end']}, Strand={orf['strand']}, "
                       f"Len={orf['length_bp']}bp, ProtLen={len(orf['protein_sequence'])}aa, Prot={orf['protein_sequence']}")
        else:
            print(f"No ORFs found [min_len={min_len_rev} aa]. (Expected)")
    except Exception as e:
        print(f"\nERROR during reverse ORF test:")
        traceback.print_exc()