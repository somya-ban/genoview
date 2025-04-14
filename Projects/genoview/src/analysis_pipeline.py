# src/analysis_pipeline.py
"""
Main Analysis Pipeline Script for GenoView.

This script orchestrates the entire sequence analysis workflow:
1. Parses an input FASTA file.
2. Finds Open Reading Frames (ORFs).
3. Detects predefined sequence motifs.
4. Identifies protein domains in predicted ORFs via InterProScan API.
5. Generates a biological summary using a Large Language Model (LLM).
6. Consolidates all results into a structured JSON output.

It can be run from the command line, taking a FASTA file as input and
optionally saving the results to a specified JSON file.
"""

import json
import os
import argparse
import datetime
import traceback # For detailed error reporting
from typing import Dict, Any, Optional # For type hinting

# --- Core Analysis Modules ---
# Import functions from sibling modules in the 'src' directory
try:
    from fasta_parser import parse_fasta_file
    from orf_finder import find_orfs_biopython
    from motif_detector import find_motifs
    from domain_finder import find_domains_interpro
except ImportError as e:
    print(f"Fatal Error: Could not import core analysis modules.\n{e}\n"
          f"Ensure all required .py files (fasta_parser, orf_finder, etc.) "
          f"are present in the 'src' directory.")
    raise # Stop execution if core components are missing

# --- LLM Integration Modules ---
# Import functions from the llm_integration sub-package
try:
    from llm_integration.prompt_builder import build_llm_prompt
    from llm_integration.summarizer import generate_huggingface_summary # Using HF summarizer
except ImportError as e:
    print(f"Fatal Error: Could not import LLM integration modules.\n{e}\n"
          f"Ensure the 'src/llm_integration' directory and its .py files exist.")
    raise # Stop execution if LLM components are missing


def run_full_analysis(fasta_filepath: str, output_json_path: Optional[str] = None, skip_llm: bool = False) -> Dict[str, Any]:
    """
    Executes the complete GenoView analysis pipeline on a single FASTA sequence.

    Handles file parsing, ORF finding, motif detection, domain searching via API,
    LLM summarization, and result aggregation.

    Args:
        fasta_filepath: Path to the input FASTA file (should contain one sequence
                        for current implementation).
        output_json_path: Optional path to save the final results as a JSON file.
                          If None, results are not saved to disk.
        skip_llm: If True, the LLM summarization step will be skipped. Useful for
                  faster runs or when API keys/quota are unavailable. Defaults to False.

    Returns:
        A dictionary containing the consolidated analysis results. Includes keys like:
        - 'sequence_id' (str)
        - 'sequence_length' (int)
        - 'analysis_timestamp' (str)
        - 'llm_summary' (str | None): The generated summary or status message.
        - 'results' (dict): Contains nested results for 'orfs', 'motifs', 'domains'.
        Returns an empty dictionary if critical errors occur early (e.g., file not found).
    """
    print(f"Starting analysis pipeline for: {fasta_filepath}")
    # Initialize the main results dictionary
    analysis_results: Dict[str, Any] = {
        "sequence_id": None,
        "sequence_length": None,
        "analysis_timestamp": datetime.datetime.utcnow().isoformat(timespec='seconds') + "Z",
        "llm_summary": None,
        "results": {
            "orfs": [],
            "motifs": [],
            "domains": []
        }
    }

    # --- 1. Parse FASTA ---
    # This step is critical, so failure here returns an empty result immediately.
    try:
        print("\nStep 1: Parsing FASTA file...")
        # Use the imported parser function
        record_iterator = parse_fasta_file(fasta_filepath)
        # Attempt to get the first sequence record
        sequence_record = next(record_iterator, None)

        # Handle case where file has no sequences
        if sequence_record is None:
            print(f"Error: No sequence record found in '{fasta_filepath}'.")
            # Return immediately as subsequent steps depend on the sequence
            return {} # Indicate critical failure

        # Check if there are more sequences in the file (current limitation)
        if next(record_iterator, None) is not None:
            print(f"Warning: Multi-FASTA file detected. GenoView currently processes only the first sequence: '{sequence_record.id}'.")

        # Store basic sequence information
        analysis_results["sequence_id"] = sequence_record.id
        analysis_results["sequence_length"] = len(sequence_record.seq)
        print(f"Processing Sequence ID: {sequence_record.id}, Length: {analysis_results['sequence_length']} bp")

    except FileNotFoundError:
        print(f"Fatal Error: Input FASTA file not found at '{fasta_filepath}'")
        return {} # Cannot proceed without input file
    except (ValueError, IOError, RuntimeError) as e:
        # Catch specific errors raised by the parser
        print(f"Fatal Error during FASTA parsing: {e}")
        return {}
    except Exception as e:
        # Catch any other unexpected errors during parsing
        print(f"Fatal Error: An unexpected error occurred during FASTA parsing:")
        traceback.print_exc()
        return {}

    # --- 2. Find ORFs ---
    # Wrap subsequent steps in try-except to allow the pipeline to continue
    # and report partial results if one step fails.
    try:
        print("\nStep 2: Finding ORFs...")
        # Use default min_protein_length=25 unless customized later
        orfs = find_orfs_biopython(sequence_record, min_protein_length=25)
        # Add unique IDs to each ORF dict *before* adding to results
        # This ID is used for linking domains back to the correct ORF
        for i, orf_dict in enumerate(orfs):
            orf_dict["orf_id"] = f"orf_{i+1}" # Generate simple ID
        analysis_results["results"]["orfs"] = orfs
        print(f"Found {len(orfs)} ORF(s) meeting criteria.")
    except Exception as e:
        print(f"Error during ORF finding:")
        traceback.print_exc()
        analysis_results["results"]["orfs"] = [{"error": "ORF finding failed."}] # Indicate failure in results

    # --- 3. Find Motifs ---
    try:
        print("\nStep 3: Finding Motifs...")
        # Uses default regex patterns defined in motif_detector.py
        motifs = find_motifs(sequence_record)
        analysis_results["results"]["motifs"] = motifs
        print(f"Found {len(motifs)} motif instance(s).")
    except Exception as e:
        print(f"Error during Motif finding:")
        traceback.print_exc()
        analysis_results["results"]["motifs"] = [{"error": "Motif detection failed."}]

    # --- 4. Prepare Proteins for Domain Search ---
    # Extract valid protein sequences from the found ORFs
    proteins_to_scan = []
    try:
        # Safely iterate over the ORFs list, checking for necessary keys
        for orf in analysis_results["results"].get("orfs", []):
            if isinstance(orf, dict) and "orf_id" in orf and isinstance(orf.get("protein_sequence"), str) and orf["protein_sequence"]:
                proteins_to_scan.append( (orf["orf_id"], orf["protein_sequence"]) )
    except Exception as e:
         print(f"Error preparing protein sequences for domain search: {e}")
         # Continue, domain search will likely be skipped

    # --- 5. Find Domains (via InterProScan API) ---
    try:
        print("\nStep 4: Finding Domains (using InterProScan API)...")
        if proteins_to_scan:
             # Call the domain finder function (handles API calls, polling, parsing)
             domains = find_domains_interpro(proteins_to_scan)
             analysis_results["results"]["domains"] = domains
             print(f"Found {len(domains)} domain match(es).")
        else:
             # Skip if no suitable proteins were found in previous steps
             print("Skipping domain search: No valid protein sequences available.")
             # Ensure domains list is empty
             analysis_results["results"]["domains"] = []
    except Exception as e:
        print(f"Error during Domain finding (InterProScan API interaction):")
        traceback.print_exc()
        analysis_results["results"]["domains"] = [{"error": "Domain search failed."}]

    # --- 6. LLM Summarization ---
    if not skip_llm:
        print("\nStep 5: Generating LLM Summary (using Hugging Face)...")
        try:
            # Build the prompt using the results gathered so far
            # Use 'mistral' format as default, suitable for the HF summarizer
            prompt = build_llm_prompt(analysis_results, format="mistral")

            if prompt and "Error:" not in prompt:
                # Call the Hugging Face summarizer function
                llm_summary = generate_huggingface_summary(prompt)

                if llm_summary:
                    print("LLM Summary generated successfully.")
                    analysis_results["llm_summary"] = llm_summary
                else:
                    # Summarizer function returned None or empty string
                    print("Failed to generate LLM summary (check summarizer logs/API status).")
                    analysis_results["llm_summary"] = "Error: Failed to generate summary via LLM."
            else:
                 # Prompt building failed
                 print(f"Skipping LLM summary: Could not build prompt ({prompt}).")
                 analysis_results["llm_summary"] = "Error: Could not build prompt for LLM."
        except Exception as e:
            # Catch any unexpected error during the summarization step
            print(f"Error during LLM Summarization step:")
            traceback.print_exc()
            analysis_results["llm_summary"] = "Error: Summarization step failed unexpectedly."
    else:
        # If skip_llm flag is True
        print("\nStep 5: Skipping LLM Summary as requested.")
        analysis_results["llm_summary"] = "Skipped by user request."

    # --- 7. Save JSON Output (Optional) ---
    if output_json_path:
        print(f"\nSaving final results to: {output_json_path}")
        try:
            # Ensure the output directory exists before trying to write the file
            output_dir = os.path.dirname(output_json_path)
            # Create directory only if a path was specified (output_dir is not empty)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                # print(f"Ensured output directory exists: {output_dir}") # Debug log

            # Write the results dictionary to the JSON file with indentation
            with open(output_json_path, 'w', encoding='utf-8') as f:
                json.dump(analysis_results, f, indent=2, ensure_ascii=False)
            print("Results saved successfully.")
        except IOError as e:
            # Handle file system errors (e.g., permissions)
            print(f"Error: Could not write JSON results to '{output_json_path}': {e}")
        except TypeError as e:
             # Handle errors if data structure is not JSON serializable
             print(f"Error: Could not serialize results to JSON: {e}")
             # Avoid printing potentially huge results dict; log snippet if needed
             # print("Problematic data structure keys:", analysis_results.keys())
        except Exception as e:
             print(f"Error: An unexpected error occurred while saving results:")
             traceback.print_exc()

    # --- Pipeline Finish ---
    print("\nAnalysis pipeline finished.")
    return analysis_results


# --- Command-Line Interface Execution ---
if __name__ == "__main__":
    # Set up argument parser for command-line usage
    parser = argparse.ArgumentParser(
        description="Run GenoView analysis pipeline on a FASTA sequence.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Show default values in help
    )
    # Required argument: Input FASTA file path
    parser.add_argument(
        "fasta_file",
        help="Path to the input FASTA file (containing a single sequence)."
    )
    # Optional argument: Output JSON file path
    parser.add_argument(
        "-o", "--output",
        help="Path to save the output JSON results. If omitted, results are printed to console.",
        default=None # Default is not to save to file
    )
    # Optional flag: Skip the LLM step
    parser.add_argument(
        "--skip-llm",
        action="store_true", # Makes it a flag, presence means True
        help="Skip the LLM summarization step (saves time/API calls)."
    )

    # Parse the arguments provided by the user
    args = parser.parse_args()

    # Basic input validation: Check if FASTA file exists
    if not os.path.exists(args.fasta_file):
         print(f"Fatal Error: Input FASTA file not found: {args.fasta_file}")
         # Exit script if input file is missing
         exit(1) # Use non-zero exit code to indicate error

    # --- Execute the Main Pipeline ---
    # Pass arguments to the core function
    final_results = run_full_analysis(args.fasta_file, args.output, args.skip_llm)

    # --- Post-Execution Summary (if not saving to file) ---
    if final_results and not args.output:
        # Print a concise summary to the console if results were generated
        # and no output file was specified.
        print("\n--- Analysis Summary (Output to Console) ---")
        print(f"Sequence ID: {final_results.get('sequence_id', 'N/A')}")
        print(f"Length: {final_results.get('sequence_length', 'N/A')} bp")
        # Safely access nested results
        results_summary = final_results.get('results', {})
        print(f"ORFs found: {len(results_summary.get('orfs', []))}")
        print(f"Motifs found: {len(results_summary.get('motifs', []))}")
        print(f"Domains found: {len(results_summary.get('domains', []))}")
        print(f"LLM Summary: {final_results.get('llm_summary', 'N/A')}")
        print("-" * 40)
        print("To see full details, re-run with the -o <filename.json> option.")

    elif not final_results:
         # If the pipeline returned an empty dictionary (due to critical early error)
         print("\nAnalysis pipeline failed to produce results.")
         exit(1) # Exit with error code
    # else: Results were saved to file, confirmation message printed within run_full_analysis

    # Exit script normally
    exit(0)