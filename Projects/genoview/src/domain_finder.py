import requests
import time
import json
import os
import traceback
from dotenv import load_dotenv # To load environment variables like email

# --- Configuration ---
# Load environment variables from a .env file if it exists
load_dotenv()

# Get required email from environment variables
# IMPORTANT: Create a file named '.env' in your project root directory
#            and add a line like: EBI_EMAIL=your_email@example.com
#            Add .env to your .gitignore file! NEVER commit your email directly.
USER_EMAIL = os.getenv("EBI_EMAIL")

if not USER_EMAIL:
    print("Warning: EBI_EMAIL environment variable not set.")
    print("InterProScan submissions require an email address.")
    # Decide how to handle: raise error, use placeholder (might fail), etc.
    # For now, let's allow proceeding but it will likely fail API-side.
    USER_EMAIL = "email_not_set@example.com" # Placeholder

# --- InterProScan API Endpoints (Verify these from current EBI documentation!) ---
# These might change. Double-check the official InterProScan REST API docs.
INTERPROSCAN_SUBMIT_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
INTERPROSCAN_STATUS_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/" # Append job_id
INTERPROSCAN_RESULT_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/" # Append job_id/xml or /json etc.

# --- Helper Function to Submit and Poll ---

def submit_interproscan_job(protein_sequence: str, sequence_id: str) -> str | None:
    """Submits a single protein sequence to InterProScan and returns the job ID."""
    payload = {
        'email': USER_EMAIL,
        'title': f'genoview_{sequence_id}', # Optional title for the job
        'goterms': 'false', # Don't need GO terms for now
        'pathways': 'false', # Don't need pathway info for now
        'stype': 'p', # Sequence type protein
        'sequence': protein_sequence
    }
    headers = {'Accept': 'text/plain'} # API typically returns Job ID as plain text on success

    try:
        response = requests.post(INTERPROSCAN_SUBMIT_URL, data=payload, headers=headers, timeout=60) # Added timeout
        response.raise_for_status() # Raise HTTPError for bad responses (4xx or 5xx)
        job_id = response.text.strip()
        print(f"  > Submitted {sequence_id}, Job ID: {job_id}")
        return job_id
    except requests.exceptions.RequestException as e:
        print(f"  > Error submitting InterProScan job for {sequence_id}: {e}")
        return None
    except Exception as e:
        print(f"  > An unexpected error occurred during job submission for {sequence_id}: {e}")
        return None

def check_job_status(job_id: str) -> str | None:
    """Checks the status of an InterProScan job."""
    try:
        url = INTERPROSCAN_STATUS_URL + job_id
        headers = {'Accept': 'text/plain'}
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
        status = response.text.strip()
        return status
    except requests.exceptions.RequestException as e:
        print(f"  > Error checking status for Job ID {job_id}: {e}")
        return "ERROR" # Return custom status for error
    except Exception as e:
        print(f"  > An unexpected error occurred during status check for {job_id}: {e}")
        return "ERROR"

def get_job_results(job_id: str) -> dict | None:
    """Retrieves the results of a finished InterProScan job in JSON format."""
    try:
        # Adjust result type if needed (e.g., /xml, /tsv) based on API docs
        url = INTERPROSCAN_RESULT_URL + job_id + "/json"
        headers = {'Accept': 'application/json'}
        response = requests.get(url, headers=headers, timeout=120) # Longer timeout for results
        response.raise_for_status()
        results = response.json()
        return results
    except json.JSONDecodeError:
        print(f"  > Error: Could not decode JSON results for Job ID {job_id}. Response was:")
        print(response.text[:500] + "...") # Print first 500 chars of response
        return None
    except requests.exceptions.RequestException as e:
        print(f"  > Error retrieving results for Job ID {job_id}: {e}")
        return None
    except Exception as e:
        print(f"  > An unexpected error occurred during result retrieval for {job_id}: {e}")
        return None

# --- Main Domain Finder Function ---

def find_domains_interpro(protein_sequences: list[tuple[str, str]], poll_interval: int = 15, max_wait_minutes: int = 30) -> list[dict]:
    """
    Finds protein domains for multiple protein sequences using the InterProScan REST API.

    Args:
        protein_sequences (list[tuple[str, str]]):
            A list of tuples, where each tuple contains:
            (orf_id: str, protein_sequence: str)
            These typically come from the output of the orf_finder.
        poll_interval (int): Seconds to wait between status checks.
        max_wait_minutes (int): Maximum time to wait for all jobs to complete.

    Returns:
        list[dict]: List of found domains, parsed from InterProScan results.
                    Keys might include: 'orf_id', 'source_db', 'accession',
                    'description', 'start_aa', 'end_aa', 'evalue', 'interpro_id', 'interpro_desc'.
    """
    all_domain_results = []
    job_ids = {} # Dictionary to track {orf_id: job_id}
    completed_jobs = {} # Dictionary for {job_id: results}
    start_time = time.time()

    print(f"\nStarting InterProScan domain finding for {len(protein_sequences)} protein sequences...")
    print(f"Using email: {USER_EMAIL}") # Confirm email being used

    if not protein_sequences:
        print("No protein sequences provided.")
        return []

    # 1. Submit all jobs
    for orf_id, protein_seq in protein_sequences:
        if not protein_seq or len(protein_seq) < 10: # Basic check for validity/length
             print(f"  > Skipping short/empty sequence: {orf_id}")
             continue
        job_id = submit_interproscan_job(protein_seq, orf_id)
        if job_id:
            job_ids[orf_id] = job_id
        time.sleep(1) # Small delay between submissions

    if not job_ids:
        print("No jobs were submitted successfully.")
        return []

    print(f"\n{len(job_ids)} jobs submitted. Polling for results (Max wait: {max_wait_minutes} mins)...")

    # 2. Poll for status and retrieve results
    submitted_job_ids = list(job_ids.values())
    while submitted_job_ids:
        current_job_id = submitted_job_ids.pop(0) # Get next job ID to check

        status = check_job_status(current_job_id)

        if status == "FINISHED":
            print(f"  > Job {current_job_id} FINISHED. Retrieving results...")
            results = get_job_results(current_job_id)
            if results:
                completed_jobs[current_job_id] = results
            else:
                # Consider job failed if results can't be retrieved
                 print(f"  > Failed to retrieve results for finished job {current_job_id}.")
                 completed_jobs[current_job_id] = None # Mark as failed/no results
        elif status == "RUNNING" or status == "QUEUED":
            # Put back in list to check again later
            submitted_job_ids.append(current_job_id)
            print(f"  > Job {current_job_id} status: {status}") # Optional: reduce verbosity here
        elif status == "ERROR" or status == "FAILURE" or status == "NOT_FOUND":
            print(f"  > Job {current_job_id} status: {status}. Will not check again.")
            completed_jobs[current_job_id] = None # Mark as failed
        else:
             print(f"  > Job {current_job_id} status: {status}. Treating as error.")
             completed_jobs[current_job_id] = None # Mark as failed


        # Check elapsed time
        elapsed_time = time.time() - start_time
        if elapsed_time > max_wait_minutes * 60:
            print(f"Error: Reached maximum wait time ({max_wait_minutes} minutes). Aborting.")
            # Mark remaining jobs as failed?
            for remaining_job_id in submitted_job_ids:
                 completed_jobs[remaining_job_id] = None
            break # Exit the polling loop

        # Wait before next check only if list is repopulated (avoid sleep if checking last item)
        if submitted_job_ids:
            print(f"  ...waiting {poll_interval}s before next check ({len(submitted_job_ids)} jobs remaining)...")
            time.sleep(poll_interval)

    print("\nPolling finished. Parsing results...")

    # 3. Parse results from completed jobs
    for orf_id, job_id in job_ids.items():
        raw_results = completed_jobs.get(job_id)
        if raw_results and 'results' in raw_results and raw_results['results']:
            # Structure of InterProScan JSON can be nested. Inspect a real example!
            # Typically results[0]['matches'] contains the list of domain hits.
            matches = raw_results['results'][0].get('matches', [])
            for match in matches:
                signature = match.get('signature', {})
                accession = signature.get('accession')
                name = signature.get('name', 'N/A')
                description = signature.get('description', 'N/A')
                db = signature.get('signatureLibrary', 'N/A').upper() # e.g., PFAM, SMART
                entry = signature.get('entry', {})
                interpro_id = entry.get('accession') if entry else None
                interpro_desc = entry.get('name') if entry else None

                # Locations are often in a list
                locations = match.get('locations', [])
                if locations:
                    # Assuming only one location per match for simplicity here
                    # Real results might have fragmented domains
                    loc = locations[0]
                    start_aa = loc.get('start')
                    end_aa = loc.get('end')
                    evalue = loc.get('evalue', 'N/A') # E-value might be here

                    # Add the parsed domain info to our list
                    all_domain_results.append({
                        "orf_id": orf_id,
                        "source_db": db,
                        "accession": accession,
                        "description": name or description, # Prefer name if available
                        "start_aa": start_aa,
                        "end_aa": end_aa,
                        "evalue": evalue,
                        "interpro_id": interpro_id, # Integrated entry ID
                        "interpro_desc": interpro_desc # Integrated entry description
                    })
        elif raw_results is None and job_id in completed_jobs:
             print(f"  > No results parsed for {orf_id} (Job {job_id} failed or timed out).")
        # else: job might not have been submitted successfully or results were empty


    print(f"Domain finding complete. Found {len(all_domain_results)} domain matches.")
    return all_domain_results

# --- Example Usage (Placeholder - Not meant for direct execution of full logic) ---
if __name__ == "__main__":
    # This block will only run if the script is executed directly.
    # It's now just a placeholder to show conceptual usage or for adding
    # very simple, non-API-calling tests in the future if needed.

    print("--- src/domain_finder.py ---")
    print("This script contains functions for finding domains via InterProScan API.")
    print("It is intended to be imported by the main analysis pipeline.")
    print("Running this script directly does not perform a full analysis.")

    # Example of how the main function *would* be called by the pipeline
    # (Do not uncomment unless testing small, specific things without API calls)
    # example_proteins_conceptual = [
    #     ("orf1_conceptual", "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"),
    # ]
    # print("\nConceptual call to find_domains_interpro (API calls disabled in this block):")
    # try:
    #      # Normally you would call find_domains_interpro here if you wanted
    #      # to test *without* relying on the pipeline script, but ensure you handle
    #      # API keys and potential costs/errors appropriately.
    #      # domains = find_domains_interpro(example_proteins_conceptual)
    #      # print(json.dumps(domains, indent=2))
    #      print("Example usage: Call find_domains_interpro(proteins_list) from another script.")
    # except NameError: # If find_domains_interpro isn't defined when running directly (shouldn't happen)
    #      print("Error: find_domains_interpro function not found.")
    # except Exception as e:
    #      print(f"An error occurred: {e}")