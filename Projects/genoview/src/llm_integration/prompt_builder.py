# src/llm_integration/prompt_builder.py
"""
Constructs prompts suitable for LLM-based summarization of sequence analysis results.

This module takes the structured dictionary output from the analysis pipeline
and formats it into a clear, informative text prompt, including specific
instructions for the desired summary format and content.
"""

import json
from typing import Dict, Any, List # For type hinting

def format_results_for_prompt(analysis_data: Dict[str, Any]) -> str:
    """
    Formats the core analysis results into a structured text block for the LLM context.

    Extracts key findings about ORFs, motifs, and domains from the pipeline output
    and presents them in a human-readable format.

    Args:
        analysis_data: The dictionary containing results from the analysis pipeline.
                       Expected to have keys like 'sequence_id', 'sequence_length',
                       and a nested 'results' dictionary with 'orfs', 'motifs', 'domains'.

    Returns:
        A formatted multi-line string summarizing the key findings.
        Returns an error message string if the input data is invalid or missing
        essential components.
    """
    # --- Input Validation ---
    if not analysis_data or "results" not in analysis_data:
        return "Error: Invalid or empty analysis data provided."

    # --- Extract Key Information ---
    seq_id: str = analysis_data.get("sequence_id", "Unknown Sequence")
    seq_len: Any = analysis_data.get("sequence_length", "N/A") # Allow N/A
    results: Dict[str, Any] = analysis_data.get("results", {})
    orfs: List[Dict[str, Any]] = results.get("orfs", [])
    motifs: List[Dict[str, Any]] = results.get("motifs", [])
    domains: List[Dict[str, Any]] = results.get("domains", [])

    # --- Build the Context String Parts ---
    context_parts: List[str] = []
    context_parts.append(f"Analysis Results for Sequence: {seq_id} (Length: {seq_len} bp)")
    context_parts.append("-" * 20) # Visual separator

    # --- ORF Summary ---
    num_orfs: int = len(orfs)
    context_parts.append(f"ORFs Found: {num_orfs}")
    if orfs:
        longest_orf_len_bp: int = 0
        longest_orf_prot_len_aa: int = 0
        try:
            # Find the ORF with the longest *protein* sequence safely
            orfs_with_prot = [orf for orf in orfs if isinstance(orf.get("protein_sequence"), str)]
            if orfs_with_prot:
                longest_orf = max(orfs_with_prot, key=lambda o: len(o.get("protein_sequence", "")))
                longest_orf_len_bp = longest_orf.get("length_bp", 0)
                longest_orf_prot_len_aa = len(longest_orf.get("protein_sequence", ""))
        except Exception as e:
            # Log potential error during calculation without stopping
            print(f"Minor error calculating longest ORF: {e}")

        if longest_orf_prot_len_aa > 0:
            context_parts.append(f"  - Longest ORF translates to {longest_orf_prot_len_aa} amino acids ({longest_orf_len_bp} bp).")
        # Future enhancement: Could add info about ORFs above a certain length threshold.

    # --- Motif Summary ---
    num_motifs: int = len(motifs)
    context_parts.append(f"Motifs Found: {num_motifs}")
    if motifs:
        # Extract unique motif IDs found, limit to first 5 for brevity in prompt
        unique_motif_ids: List[str] = sorted(list(set(m.get("motif_id", "N/A") for m in motifs)))
        display_motifs: str = ", ".join(unique_motif_ids[:5])
        if len(unique_motif_ids) > 5:
            display_motifs += "..."
        context_parts.append(f"  - Types detected include: {display_motifs}")

    # --- Domain Summary (Often the most informative part) ---
    num_domains: int = len(domains)
    context_parts.append(f"Protein Domains Found: {num_domains}")
    if domains:
        domain_summary_lines: List[str] = []
        # Limit the number of detailed domain entries shown in the prompt
        max_domains_in_prompt: int = 5
        displayed_count: int = 0

        for d in domains:
            # Extract details safely using .get() with defaults
            db: str = d.get("source_db", "N/A")
            acc: str = d.get("accession", "N/A")
            desc: str = d.get("description", "N/A")
            orf_id: str = d.get("orf_id", "N/A") # Link back to the parent ORF
            interpro: str = f" (InterPro: {d.get('interpro_id', 'N/A')})" if d.get('interpro_id') else ""

            if displayed_count < max_domains_in_prompt:
                # Format the domain line clearly
                domain_summary_lines.append(f"  - {db}:{acc} ({desc}){interpro} [found in {orf_id}]")
                displayed_count += 1
            elif displayed_count == max_domains_in_prompt:
                # Add ellipsis once if more domains exist
                domain_summary_lines.append("  - ... (additional domains found)")
                displayed_count += 1 # Prevent adding ellipsis repeatedly

        context_parts.extend(domain_summary_lines)

    # --- Combine all parts into a single string ---
    return "\n".join(context_parts)


def build_llm_prompt(analysis_data: Dict[str, Any], format: str = "mistral") -> str:
    """
    Builds the full prompt string for the LLM, combining context and instructions.

    Supports different formatting conventions (e.g., 'basic', 'mistral') based
    on the target LLM's requirements or preferences.

    Args:
        analysis_data: The dictionary containing results from the analysis pipeline.
        format: The desired prompt format. Currently supports 'basic' (generic)
                and 'mistral' (using [INST] tags). Defaults to 'mistral'.

    Returns:
        The complete prompt string ready to be sent to the LLM API.
        Returns a fallback prompt if context generation fails.
    """
    # Generate the context summary text block
    context_summary: str = format_results_for_prompt(analysis_data)

    # Handle potential errors during context formatting
    if "Error:" in context_summary:
        print(f"Warning: Could not format results for prompt generation: {context_summary}")
        # Provide a generic fallback prompt if context is unavailable
        return "Please provide a summary of general biological sequence analysis."

    # --- Define the Core Instruction for the LLM ---
    # Be specific about the desired output (length, focus, style)
    instruction_core: str = """Based on the provided analysis results, please generate a concise biological summary (2-3 sentences maximum).
Focus on the most significant findings, particularly any identified protein domains and their potential functional implications.
If ORFs were found, mention if they are likely protein-coding based on domain hits. If no significant features were found, state that clearly.
Interpret the findings briefly; do not simply list the results again."""

    # --- Assemble the Final Prompt Based on Format ---
    full_prompt: str
    if format == "mistral":
        # Format specifically for Mistral Instruct models
        # Structure: <s>[INST] Instruction + Context [/INST] Desired Start of Output
        full_prompt = f"<s>[INST] {instruction_core}\n\nAnalysis Context:\n{context_summary} [/INST]Biological Summary:"
    else: # Default to a basic, generally compatible format
        full_prompt = f"""Analysis Context:
{context_summary}

Instructions:
{instruction_core}

Biological Summary:
"""
    # The final "Biological Summary:" acts as a prompt for the model to start generating

    return full_prompt

# --- Example Usage Block ---
if __name__ == "__main__":
    print("--- Testing prompt_builder.py ---")

    # Define sample analysis data dictionaries for testing
    dummy_data_with_results = {
        "sequence_id": "GeneX_Example", "sequence_length": 1500,
        "results": {
            "orfs": [{"orf_id": "orf_1", "start": 100, "end": 1300, "strand": "+", "length_bp": 1201, "protein_sequence": "MKT..."},
                     {"orf_id": "orf_2", "start": 1400, "end": 1480, "strand": "-", "length_bp": 81, "protein_sequence": "MSHORT"}],
            "motifs": [{"motif_id": "TATA_BOX_LIKE", "start": 30, "end": 37, "strand": "+", "matched_sequence": "TATAATAA"},
                       {"motif_id": "GC_BOX_LIKE", "start": 70, "end": 75, "strand": "+", "matched_sequence": "GGGCGG"}],
            "domains": [{"orf_id": "orf_1", "source_db": "PFAM", "accession": "PF00069", "description": "Protein kinase domain", "start_aa": 50, "end_aa": 300, "evalue": "1e-60", "interpro_id": "IPR000719", "interpro_desc": "Protein kinase domain"},
                        {"orf_id": "orf_1", "source_db": "SMART", "accession": "SM00220", "description": "Ser/Thr protein kinase active site", "start_aa": 150, "end_aa": 170, "evalue": "5e-15", "interpro_id": "IPR017442", "interpro_desc": "Protein kinase, active site"}]
        }
    }
    dummy_data_no_results = {
        "sequence_id": "Short_Example", "sequence_length": 90,
        "results": {"orfs": [], "motifs": [], "domains": []}
    }

    # Test generating prompts in different formats
    print("\n--- Prompt (With Results, Mistral Format) ---")
    prompt1_mistral = build_llm_prompt(dummy_data_with_results, format="mistral")
    print(prompt1_mistral)

    print("\n--- Prompt (With Results, Basic Format) ---")
    prompt1_basic = build_llm_prompt(dummy_data_with_results, format="basic")
    print(prompt1_basic)

    print("\n--- Prompt (No Results, Mistral Format) ---")
    prompt2_mistral = build_llm_prompt(dummy_data_no_results, format="mistral")
    print(prompt2_mistral)