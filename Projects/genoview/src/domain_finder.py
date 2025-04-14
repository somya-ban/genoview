# src/domain_finder.py
"""
Finds protein domain annotations using the InterProScan REST API.

This module handles submitting protein sequences derived from predicted ORFs
to the EBI InterProScan web service, polling for job completion, retrieving,
and parsing the results to identify conserved protein domains and signatures.

Requires internet connectivity and a valid email address provided via the
EBI_EMAIL environment variable (in a .env file) for API usage.
"""

import requests
import time
import json
import os
import traceback
from dotenv import load_dotenv # To load environment variables
from typing import List, Tuple, Dict, Any, Optional # For type hinting

# --- Configuration ---
# Load environment variables from a .env file in the project root
load_dotenv()

# Get required email from environment variable for EBI API calls
# It's crucial for tracking and compliance with EBI policies.
USER_EMAIL: Optional[str] = os.getenv("EBI_EMAIL")

# Check if email is set, warn if not. API calls likely fail without it.
if not USER_EMAIL:
    print("CRITICAL WARNING: EBI_EMAIL environment variable not set in .env file.")
    print("InterProScan submissions WILL LIKELY FAIL without a valid email.")
    # Proceeding might still allow testing other parts, but API interaction is compromised.

# InterProScan v5 REST API Endpoints
# IMPORTANT: These URLs are subject to change. Always verify with the latest
#            official EBI/InterProScan web services documentation.
#            (e.g., search "InterProScan REST API Documentation")
INTERPROSCAN_SUBMIT_URL: str = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
INTERPROSCAN_STATUS_URL: str = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/" # Append job_id
INTERPROSCAN_RESULT_URL: str = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/" # Append job_id/result_type

# --- API Interaction Helper Functions ---

def submit_interproscan_job(protein_sequence: str, sequence_id: str) -> Optional[str]:
    """
    Submits a single protein sequence to the InterProScan REST API.

    Args:
        protein_sequence: The amino acid sequence string to be analyzed.
        sequence_id: A unique identifier for this sequence (e.g., derived from ORF ID).

    Returns:
        The Job ID string if submission is successful, otherwise None.
    """
    if not USER_EMAIL:
        print(f"  > Error: Cannot submit {sequence_id}, EBI_EMAIL not configured.")
        return None

    # Define payload parameters for the InterProScan job submission
    payload = {
        'email': USER_EMAIL,          # Required user email
        'title': f'genoview_{sequence_id}', # Optional job title
        'goterms': 'false',           # Exclude GO Term enrichment for faster results
        'pathways': 'false',          # Exclude pathway information
        'stype': 'p',                 # Sequence type is protein ('p')
        'sequence': protein_sequence  # The actual protein sequence
        # Add 'appl' parameter to specify databases if needed, e.g., 'appl': 'Pfam,SMART'
    }
    # Specify expected response type (Job ID is plain text)
    headers = {'Accept': 'text/plain'}
    # Set a reasonable timeout for the submission request
    submit_timeout: int = 60 # seconds

    try:
        response = requests.post(INTERPROSCAN_SUBMIT_URL, data=payload, headers=headers, timeout=submit_timeout)
        # Check for HTTP errors (e.g., 4xx client errors, 5xx server errors)
        response.raise_for_status()
        # Successful submission returns the Job ID as plain text
        job_id: str = response.text.strip()
        print(f"  > Submitted {sequence_id}, Job ID: {job_id}")
        return job_id
    except requests.exceptions.Timeout:
        print(f"  > Error: Timeout submitting InterProScan job for {sequence_id}.")
        return None
    except requests.exceptions.RequestException as e:
        print(f"  > Error submitting InterProScan job for {sequence_id}: {e}")
        # Print response body if available for more details
        if e.response is not None:
             print(f"  > Response status: {e.response.status_code}, Body: {e.response.text[:200]}...")
        return None
    except Exception as e:
        # Catch any other unexpected errors during submission
        print(f"  > An unexpected error occurred during job submission for {sequence_id}:")
        traceback.print_exc()
        return None

def check_job_status(job_id: str) -> Optional[str]:
    """
    Checks the status (e.g., RUNNING, FINISHED, ERROR) of an InterProScan job.

    Args:
        job_id: The Job ID string returned by the submission step.

    Returns:
        The status string (e.g., "RUNNING", "FINISHED", "ERROR", "FAILURE",
        "NOT_FOUND") if successful, otherwise None if the check fails.
    """
    try:
        url: str = INTERPROSCAN_STATUS_URL + job_id
        headers = {'Accept': 'text/plain'}
        status_timeout: int = 30 # seconds

        response = requests.get(url, headers=headers, timeout=status_timeout)
        response.raise_for_status()
        status: str = response.text.strip().upper() # Standardize status to uppercase
        return status
    except requests.exceptions.Timeout:
         print(f"  > Error: Timeout checking status for Job ID {job_id}.")
         return "ERROR" # Treat timeout during status check as an error state
    except requests.exceptions.RequestException as e:
        print(f"  > Error checking status for Job ID {job_id}: {e}")
        # Return specific status if possible, else generic error
        if e.response is not None and e.response.status_code == 404:
             return "NOT_FOUND"
        return "ERROR" # Treat other request errors as error state
    except Exception as e:
        print(f"  > An unexpected error occurred during status check for {job_id}:")
        traceback.print_exc()
        return "ERROR"

def get_job_results(job_id: str, result_type: str = "json") -> Optional[Dict[str, Any]]:
    """
    Retrieves the results of a finished InterProScan job.

    Args:
        job_id: The Job ID string of a job confirmed to be "FINISHED".
        result_type: The desired format ('json', 'xml', 'tsv', 'gff').
                     Defaults to 'json'.

    Returns:
        A dictionary containing the parsed results (if JSON requested and
        successful), otherwise None.
    """
    try:
        url: str = INTERPROSCAN_RESULT_URL + job_id + "/" + result_type
        # Request JSON format
        headers = {'Accept': 'application/json' if result_type == "json" else '*/*'}
        result_timeout: int = 180 # Allow more time for potentially large results

        response = requests.get(url, headers=headers, timeout=result_timeout)
        response.raise_for_status()

        # Parse based on requested type
        if result_type == "json":
            results_data: Dict[str, Any] = response.json()
            return results_data
        else:
            # For other types, could return raw text or implement specific parsers
            print(f"Warning: Result type '{result_type}' requested. Returning raw text.")
            # For simplicity, we primarily handle JSON here. Return None for others.
            # return response.text
            return None

    except json.JSONDecodeError:
        # Handle cases where the response wasn't valid JSON
        print(f"  > Error: Could not decode JSON results for Job ID {job_id}. Response was:")
        print(f"  > {response.text[:500]}...") # Show beginning of invalid response
        return None
    except requests.exceptions.Timeout:
         print(f"  > Error: Timeout retrieving results for Job ID {job_id}.")
         return None
    except requests.exceptions.RequestException as e:
        print(f"  > Error retrieving results for Job ID {job_id}: {e}")
        if e.response is not None:
             print(f"  > Response status: {e.response.status_code}, Body: {e.response.text[:200]}...")
        return None
    except Exception as e:
        print(f"  > An unexpected error occurred during result retrieval for {job_id}:")
        traceback.print_exc()
        return None

def parse_interproscan_json(raw_results: Dict[str, Any], orf_id: str) -> List[Dict[str, Any]]:
    """
    Parses the domain/signature matches from the InterProScan JSON results structure.

    Args:
        raw_results: The dictionary obtained from get_job_results(..., result_type="json").
        orf_id: The original ORF ID this result corresponds to, for linking.

    Returns:
        A list of dictionaries, each representing a found domain/signature match.
        Returns an empty list if no matches are found or parsing fails.
    """
    parsed_domains: List[Dict[str, Any]] = []
    if not raw_results or 'results' not in raw_results or not raw_results['results']:
        # Expected structure: results is a list, usually with one element
        return []

    try:
        # InterProScan JSON structure can be complex. This assumes common patterns.
        # The main results are often in the first element of the 'results' list.
        matches = raw_results['results'][0].get('matches', [])

        for match in matches:
            signature = match.get('signature', {})
            if not signature: continue # Skip if no signature info

            accession = signature.get('accession')
            name = signature.get('name') # Often more descriptive than 'description'
            description = signature.get('description') # Can be None
            db = signature.get('signatureLibrary', 'N/A').upper() # e.g., PFAM, SMART, CDD
            entry = signature.get('entry', {}) # Contains integrated InterPro entry info
            interpro_id = entry.get('accession') if entry else None
            interpro_desc = entry.get('name') if entry else None

            # Locations contain coordinates and scores/e-values
            locations = match.get('locations', [])
            if locations:
                # A single signature match can have multiple fragments/locations
                # For simplicity, we process the first location here.
                # A more robust parser might handle multiple fragments per match.
                loc = locations[0]
                start_aa = loc.get('start')
                end_aa = loc.get('end')
                # E-value location can vary, check common places
                evalue = loc.get('evalue', match.get('score', 'N/A')) # Use match score as fallback?

                # Ensure essential fields are present
                if accession and start_aa is not None and end_aa is not None:
                    parsed_domains.append({
                        "orf_id": orf_id,
                        "source_db": db,
                        "accession": accession,
                        # Prefer 'name' for description, fallback to 'description' field
                        "description": name or description or "N/A",
                        "start_aa": int(start_aa),
                        "end_aa": int(end_aa),
                        "evalue": str(evalue), # Store E-value as string for flexibility
                        "interpro_id": interpro_id,
                        "interpro_desc": interpro_desc
                    })
            # else: Match had no location info, skip?

    except (KeyError, IndexError, TypeError, ValueError) as e:
        print(f"  > Error parsing InterProScan JSON structure for ORF {orf_id}: {e}")
        # Consider logging the problematic raw_results snippet here for debugging
        # print(f"Problematic JSON snippet: {str(raw_results)[:1000]}")
    except Exception as e:
        print(f"  > An unexpected error occurred during JSON parsing for {orf_id}:")
        traceback.print_exc()

    return parsed_domains


# --- Main Domain Finder Function ---

def find_domains_interpro(protein_sequences: List[Tuple[str, str]], poll_interval: int = 20, max_wait_minutes: int = 30) -> List[Dict[str, Any]]:
    """
    Finds protein domains for multiple protein sequences using the InterProScan REST API.

    Manages job submission, status polling, result retrieval, and parsing for a
    list of protein sequences.

    Args:
        protein_sequences: A list of tuples, where each tuple contains:
                           (orf_id: str, protein_sequence: str).
                           These typically come from the output of the orf_finder.
        poll_interval: Seconds to wait between checking job statuses. Avoid
                       setting too low to respect EBI servers. Defaults to 20s.
        max_wait_minutes: Maximum total time (in minutes) to wait for all jobs
                          to complete before aborting. Defaults to 30 mins.

    Returns:
        A list of dictionaries, each representing a found domain/signature,
        parsed from successful InterProScan results. Returns an empty list if
        no domains are found, no jobs complete successfully, or errors occur.
    """
    all_parsed_domains: List[Dict[str, Any]] = []
    # Dictionary to map submitted orf_id to the returned job_id
    submitted_jobs: Dict[str, str] = {}
    # Dictionary to store results or failure status for each job_id
    job_outcomes: Dict[str, Optional[Dict[str, Any]]] = {}
    start_time: float = time.time()

    # --- Input Validation ---
    if not protein_sequences:
        print("Info: No protein sequences provided to domain finder.")
        return []
    if not USER_EMAIL:
         print("Error: Cannot run domain search, EBI_EMAIL is not configured.")
         return []

    print(f"\nStarting InterProScan domain finding for {len(protein_sequences)} protein sequences...")
    print(f"Using email: {USER_EMAIL}") # Confirm email being used

    # --- 1. Submit All Jobs ---
    print("Submitting jobs...")
    for orf_id, protein_seq in protein_sequences:
        # Basic validation: Skip obviously invalid/short sequences
        if not protein_seq or len(protein_seq) < 5: # InterProScan might have min length
             print(f"  > Skipping short/empty sequence: {orf_id}")
             continue

        job_id = submit_interproscan_job(protein_seq, orf_id)
        if job_id:
            submitted_jobs[orf_id] = job_id
            job_outcomes[job_id] = None # Mark as submitted, outcome pending
        else:
            # Submission failed, no need to track this orf_id further
            print(f"  > Failed to submit job for {orf_id}.")

        # Add a small delay between submissions to be polite to the API
        time.sleep(1)

    if not submitted_jobs:
        print("No jobs were submitted successfully.")
        return []

    # --- 2. Poll for Status and Retrieve Results ---
    # List of job IDs we are still waiting for
    pending_job_ids: List[str] = list(submitted_jobs.values())
    print(f"\n{len(pending_job_ids)} jobs submitted. Polling for results (Interval: {poll_interval}s, Max wait: {max_wait_minutes} mins)...")

    while pending_job_ids:
        # Check for timeout
        elapsed_time_secs: float = time.time() - start_time
        if elapsed_time_secs > max_wait_minutes * 60:
            print(f"Error: Reached maximum wait time ({max_wait_minutes} minutes). Aborting polling.")
            # Mark any remaining pending jobs as timed out
            for job_id in pending_job_ids:
                 job_outcomes[job_id] = {"error": "Timeout"} # Indicate timeout
            break # Exit the polling loop

        current_job_id: str = pending_job_ids.pop(0) # Get the next job to check
        status = check_job_status(current_job_id)

        # Process based on status
        if status == "FINISHED":
            print(f"  > Job {current_job_id} FINISHED. Retrieving results...")
            results = get_job_results(current_job_id, result_type="json")
            if results:
                job_outcomes[current_job_id] = results # Store successful results
            else:
                 print(f"  > Failed to retrieve results for finished job {current_job_id}.")
                 job_outcomes[current_job_id] = {"error": "ResultRetrievalFailed"}
            # Job is complete, do not requeue

        elif status == "RUNNING" or status == "QUEUED":
            # Job is still pending, put it back at the end of the list to check later
            pending_job_ids.append(current_job_id)
            # Optional: Reduce print frequency for running jobs to avoid log spam
            # if time.time() % (poll_interval * 5) < poll_interval: # Print status every 5 polls
            #    print(f"  > Job {current_job_id} status: {status}...")

        elif status in ["ERROR", "FAILURE", "NOT_FOUND"]:
            # Job failed permanently or was not found
            print(f"  > Job {current_job_id} status: {status}. Marking as failed.")
            job_outcomes[current_job_id] = {"error": status}
            # Job is complete (failed), do not requeue

        else: # Handle unexpected status strings
             print(f"  > Job {current_job_id} has unexpected status: '{status}'. Treating as error.")
             job_outcomes[current_job_id] = {"error": f"UnexpectedStatus: {status}"}
             # Job is complete (failed), do not requeue

        # Only sleep if there are still pending jobs to check in this cycle
        if pending_job_ids:
            time.sleep(poll_interval)
        # Add a small final status print when moving to parsing
        elif not pending_job_ids:
             print("  > All jobs processed or failed.")


    print("\nPolling finished. Parsing results...")

    # --- 3. Parse Results from Completed/Retrieved Jobs ---
    for orf_id, job_id in submitted_jobs.items():
        outcome = job_outcomes.get(job_id)
        if isinstance(outcome, dict) and "error" not in outcome:
            # Successfully retrieved results dictionary
            parsed = parse_interproscan_json(outcome, orf_id)
            if parsed:
                 all_parsed_domains.extend(parsed)
            # else: print(f"  > No parsable matches found for {orf_id} (Job {job_id}).")
        elif isinstance(outcome, dict) and "error" in outcome:
            # Job failed or timed out
            print(f"  > Skipping parsing for failed/timed-out job {job_id} (ORF: {orf_id}), Error: {outcome['error']}")
        # else: Job outcome was None (shouldn't happen if logic above is correct)

    print(f"\nDomain finding complete. Found {len(all_parsed_domains)} domain matches across all sequences.")
    return all_parsed_domains


# --- Placeholder Example Usage Block ---
if __name__ == "__main__":
    # This block is primarily for documentation or very basic, non-API tests.
    # The main functionality is tested via the analysis pipeline or direct function calls.
    print("--- src/domain_finder.py ---")
    print("Contains functions for finding domains via InterProScan API.")
    print("Intended for import by the analysis pipeline.")
    print("Ensure EBI_EMAIL is set in .env for API calls.")
    # Example conceptual call:
    # test_proteins = [("test_orf", "MAGW...")]
    # domains = find_domains_interpro(test_proteins)
    # print(domains)