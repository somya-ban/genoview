from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data import CodonTable
import traceback # Keep for error handling in main block

def find_orfs_biopython(sequence_record: SeqRecord, min_protein_length: int = 25, table: int = 1) -> list[dict]:
    """
    Finds Open Reading Frames (ORFs) in a DNA sequence using Biopython.

    Identifies potential ORFs based on standard start ('M') and stop ('*') codons
    in all 6 reading frames.

    Args:
        sequence_record (SeqRecord): A Biopython SeqRecord object containing the DNA sequence.
        min_protein_length (int): Minimum length of the translated protein sequence
                                   (excluding stop codon) to be considered an ORF. Defaults to 25.
        table (int): NCBI Genetic Code table number to use for translation. Defaults to 1 (Standard Code).

    Returns:
        list[dict]: A list of dictionaries, where each dictionary represents an ORF
                    with keys: 'start', 'end', 'strand', 'length_bp',
                    'protein_sequence'. Returns an empty list if no ORFs are found or on error.
                    Coordinates are 1-based.
    """
    orf_results = []
    seq = sequence_record.seq
    seq_len = len(seq)

    try:
        codon_table = CodonTable.ambiguous_dna_by_id[table]
    except KeyError:
        print(f"Warning: Codon table {table} not found. Using standard table 1.")
        codon_table = CodonTable.ambiguous_dna_by_id[1]

    # Process both forward (+) and reverse (-) strands
    for strand, nuc_seq in [(+1, seq), (-1, seq.reverse_complement())]:
        strand_char = '+' if strand == 1 else '-'

        # Process the three reading frames (0, 1, 2) for this strand
        for frame in range(3):
            frame_seq = nuc_seq[frame:]
            protein_seq = ""
            dna_indices = [] # Store DNA start index (0-based) for each codon in this frame

            # Translate codon by codon
            for i in range(0, len(frame_seq) - (len(frame_seq) % 3), 3):
                codon = frame_seq[i:i+3]
                # Store the 0-based start index of this codon relative to the nuc_seq (fwd or revcomp)
                dna_indices.append(frame + i)
                try:
                    amino_acid = codon.translate(table=table, stop_symbol="*")
                    protein_seq += amino_acid
                except CodonTable.TranslationError:
                    protein_seq += "?" # Represent unknown translation

            if not protein_seq: continue # Skip if no protein sequence generated

            # Scan for M...* within the translated protein sequence
            prot_len = len(protein_seq)
            current_orf_start_prot = -1 # 0-based index in protein_seq where current ORF started

            for aa_index in range(prot_len):
                aa = protein_seq[aa_index]

                # If we find a start codon ('M') and are not already in an ORF
                if aa == 'M' and current_orf_start_prot == -1:
                    current_orf_start_prot = aa_index

                # If we find a stop codon ('*') and are currently tracking an ORF
                elif aa == '*' and current_orf_start_prot != -1:
                    current_orf_end_prot = aa_index
                    prot_orf_len = current_orf_end_prot - current_orf_start_prot

                    # Check if the protein ORF meets the minimum length criteria
                    if prot_orf_len >= min_protein_length:
                        # Calculate DNA coordinates (1-based) using stored indices
                        # DNA start coord (1-based) = (0-based index of start codon on nuc_seq) + 1
                        dna_orf_start_on_strand = dna_indices[current_orf_start_prot] + 1
                        # DNA end coord (1-based) = (0-based index of stop codon on nuc_seq) + 3 (length of codon)
                        dna_orf_end_on_strand = dna_indices[current_orf_end_prot] + 3

                        orf_protein = protein_seq[current_orf_start_prot:current_orf_end_prot]
                        # Calculate DNA length based on the calculated 1-based coords
                        orf_dna_len = dna_orf_end_on_strand - dna_orf_start_on_strand + 1

                        # Map back to forward strand if necessary
                        if strand == -1:
                            final_dna_start = seq_len - dna_orf_end_on_strand + 1
                            final_dna_end = seq_len - dna_orf_start_on_strand + 1
                        else:
                            final_dna_start = dna_orf_start_on_strand
                            final_dna_end = dna_orf_end_on_strand

                        orf_results.append({
                            "start": final_dna_start,
                            "end": final_dna_end,
                            "strand": strand_char,
                            "length_bp": orf_dna_len,
                            "protein_sequence": orf_protein
                        })

                    # Reset ORF tracking after finding a stop (regardless of length condition)
                    current_orf_start_prot = -1

    # Sort results by start position (optional but nice)
    orf_results.sort(key=lambda x: x['start'])
    return orf_results

# --- Example Usage (for testing this script directly) ---
if __name__ == "__main__":
    from Bio.Seq import Seq
    # Example DNA sequence with a clear ORF
    # Forward strand: Frame 0 (starts index 0): G K P F M K K K K K G P F * (Prot len 9 after M)
    # DNA Coords: Start codon ATG at 13 (1-based), Stop codon TAG at 43 (1-based)
    # ORF DNA: ATG...TTT (len 33 bp) -> End coord should be Start + len - 1 = 13+33-1=45?
    # Let's re-check: M @ prot idx 4, * @ prot idx 13.
    # DNA idx M: dna_indices[4] = 0 + (4*3) = 12. DNA coord = 12+1=13
    # DNA idx *: dna_indices[13] = 0 + (13*3) = 39. DNA coord end = 39+3=42.
    # ORF DNA len = 42 - 13 + 1 = 30 bp. Protein = 9 aa. Seems short?
    # Ah, the codons are GGG AAA CCC TTT / ATG AAA AAA AAA AAA / AAA GGG CCC TTT / TAG GGC CC
    # Frame 0 starts GGG.
    # Protein: G K P F / M K K K K / K G P F / * G P -> Stop is at index 13. Yes.
    # DNA indices: 0, 3, 6, 9, / 12, 15, 18, 21, 24, / 27, 30, 33, 36 / 39, 42, 45
    # Start M is protein index 4, DNA index 12 -> Start coord 13.
    # Stop * is protein index 13, DNA index 39 -> End coord 39+3 = 42.
    # Length bp = 42 - 13 + 1 = 30. Correct.
    test_dna = "GGGAAACCCTTTATGAAAAAAAAAAAAAAAGGGCCCTTTTAGGGCCC" # Len 47
    record_id = "test_seq_orf"
    test_record = SeqRecord(Seq(test_dna), id=record_id, description="Test sequence for ORF finding")

    print(f"Finding ORFs in sequence: {record_id} (Length: {len(test_record.seq)} bp)")
    print(f"Sequence: {test_record.seq}")

    # --- Set minimum protein length for this test ---
    min_prot_len = 9 # Minimum 9 amino acids protein length (matches the expected ORF)

    try:
        # Use the Biopython-based function
        found_orfs = find_orfs_biopython(test_record, min_protein_length=min_prot_len)

        if found_orfs:
            print(f"\nFound {len(found_orfs)} ORF(s) with minimum protein length {min_prot_len} aa:")
            for i, orf in enumerate(found_orfs):
                print(f"--- ORF {i+1} ---")
                print(f"  Start: {orf['start']}")
                print(f"  End: {orf['end']}")
                print(f"  Strand: {orf['strand']}")
                print(f"  Length: {orf['length_bp']} bp")
                print(f"  Protein Length: {len(orf['protein_sequence'])} aa")
                print(f"  Protein: {orf['protein_sequence']}")
        else:
            print(f"\nNo ORFs found with minimum protein length {min_prot_len} aa.")

    except Exception as e:
        print(f"\nAn unexpected error occurred during testing:")
        traceback.print_exc() # Print detailed traceback for debugging

    # --- Example with reverse complement ---
    print("\n" + "="*30 + "\n")
    # RevComp: CCCGGGCCTTTATGAAAAAAAAAATGTAGCCCGGG (Len 36)
    # Frame 1 (starts index 1): R A F M K K K K C S P (Prot len 8 after M, No stop)
    test_dna_rev = "CCCGGGCTACATTTTTTTTTTTCATAAAGGCCCGGG" # Len 36
    record_id_rev = "test_seq_orf_rev"
    test_record_rev = SeqRecord(Seq(test_dna_rev), id=record_id_rev, description="Test sequence for reverse ORF finding")

    print(f"Finding ORFs in sequence: {record_id_rev} (Length: {len(test_record_rev.seq)} bp)")
    print(f"Sequence: {test_record_rev.seq}")

    # --- Set minimum protein length for this test ---
    min_prot_len_rev = 5

    try:
        # Use the Biopython-based function
        found_orfs_rev = find_orfs_biopython(test_record_rev, min_protein_length=min_prot_len_rev)
        if found_orfs_rev:
            print(f"\nFound {len(found_orfs_rev)} ORF(s) with minimum protein length {min_prot_len_rev} aa:")
            for i, orf in enumerate(found_orfs_rev):
                print(f"--- ORF {i+1} ---")
                print(f"  Start: {orf['start']}")
                print(f"  End: {orf['end']}")
                print(f"  Strand: {orf['strand']}")
                print(f"  Length: {orf['length_bp']} bp")
                print(f"  Protein Length: {len(orf['protein_sequence'])} aa")
                print(f"  Protein: {orf['protein_sequence']}")
        else:
            print(f"\nNo ORFs found with minimum protein length {min_prot_len_rev} aa.") # Expected for this sequence
    except Exception as e:
        print(f"\nAn unexpected error occurred during testing:")
        traceback.print_exc() # Print detailed traceback for debugging