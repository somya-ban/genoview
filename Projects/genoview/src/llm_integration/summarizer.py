# src/llm_integration/summarizer.py
"""
Handles interaction with the LLM API to generate summaries.

This module takes a formatted prompt, sends it to the configured LLM API endpoint
(currently set up for the Hugging Face Inference API), and returns the
generated summary text. Requires a Hugging Face API token set in the .env file.
"""

import os
import requests # Using requests for Hugging Face API
from dotenv import load_dotenv
import traceback
import time # For potential retries or delays
from typing import Optional, Dict, Any # For type hinting

# Load environment variables (.env file in project root)
load_dotenv()

# --- Configuration ---
# Retrieve Hugging Face API Token from environment variables
API_TOKEN: Optional[str] = os.getenv("HF_API_TOKEN")

# Default Hugging Face Model Endpoint URL
# Mistral-7B-Instruct is generally a strong open model for this task.
# Verify model availability and endpoint on huggingface.co/inference-api
DEFAULT_HF_API_URL: str = "https://api-inference.huggingface.co/models/mistralai/Mistral-7B-Instruct-v0.2"

# Authentication Headers for Hugging Face API
if API_TOKEN:
    HEADERS: Dict[str, str] = {"Authorization": f"Bearer {API_TOKEN}"}
else:
    HEADERS = {} # Will likely cause authentication errors if token is missing
    print("CRITICAL WARNING: HF_API_TOKEN environment variable not set in .env file.")
    print("Hugging Face API calls WILL LIKELY FAIL without authentication.")


# --- Summarization Function ---
def generate_huggingface_summary(prompt: str,
                                 model_url: str = DEFAULT_HF_API_URL,
                                 max_new_tokens: int = 150,
                                 temperature: float = 0.6,
                                 wait_for_model: bool = True,
                                 timeout: int = 120) -> Optional[str]:
    """
    Generates a summary using the Hugging Face Inference API for text generation.

    Args:
        prompt: The complete input prompt string for the LLM.
        model_url: The Hugging Face Inference API URL for the desired model.
                   Defaults to `DEFAULT_HF_API_URL`.
        max_new_tokens: Max number of tokens for the *generated* summary.
        temperature: Sampling temperature (must be > 0 for HF API). Controls randomness.
                     Lower values (~0.1-0.5) are more focused, higher values (~0.7-1.0)
                     are more creative/random. Defaults to 0.6.
        wait_for_model: If True (default), the request will wait if the model
                        is currently loading on the HF infrastructure (can take
                        minutes on free tier). If False, returns an error immediately
                        if the model is not ready.
        timeout: Request timeout in seconds. Defaults to 120.

    Returns:
        The generated summary text string if successful, otherwise None.
    """
    if not API_TOKEN:
        print("Error: Hugging Face API Token not set. Cannot generate summary.")
        return None
    if not HEADERS: # Double check based on token presence
         print("Error: Authentication headers not set. Cannot generate summary.")
         return None

    # Display which model is being called
    model_name = model_url.split('/')[-1] if '/' in model_url else model_url
    print(f"\n--- Sending prompt to Hugging Face model: {model_name} ---")
    # print(f"Prompt snippet: {prompt[:200]}...") # Uncomment for debugging

    # Construct the payload according to HF Inference API documentation
    # for text-generation tasks.
    payload: Dict[str, Any] = {
        "inputs": prompt,
        "parameters": {
            # Max tokens the *model* should generate, not total length
            "max_new_tokens": max_new_tokens,
            # Ensure temperature is slightly above zero if 0.0 is passed
            "temperature": max(temperature, 0.01),
            # Only return the generated text, not the prompt input
            "return_full_text": False,
            # Request only one generated sequence
            "num_return_sequences": 1,
            # Optional: Add other parameters like 'top_p', 'top_k', 'repetition_penalty'
        },
        "options": {
             # Waits for model if loading (recommended for free/serverless tier)
             "wait_for_model": wait_for_model
        }
    }

    response: Optional[requests.Response] = None # Initialize response variable
    try:
        # Send the POST request to the Hugging Face API endpoint
        response = requests.post(model_url, headers=HEADERS, json=payload, timeout=timeout)

        # Check for HTTP errors (e.g., 401 Unauthorized, 429 Rate Limit, 503 Service Unavailable)
        response.raise_for_status()

        # Parse the JSON response
        result: Any = response.json()

        # Process the result - structure can vary slightly.
        # Usually a list containing a dictionary with "generated_text".
        if isinstance(result, list) and len(result) > 0 and isinstance(result[0], dict):
            summary: str = result[0].get("generated_text", "").strip()
            if summary:
                 print("--- Hugging Face response received ---")
                 return summary
            else:
                 # Handle case where 'generated_text' key exists but is empty
                 print("Error: Hugging Face response contained empty 'generated_text'.")
                 print("Response:", result)
                 return None
        else:
            # Handle unexpected response formats
            print("Error: Hugging Face response format unexpected.")
            print("Response:", result)
            return None

    # --- Specific Error Handling ---
    except requests.exceptions.Timeout:
        print(f"Error: Request timed out calling Hugging Face API ({model_name}). "
              f"Model might be unavailable or analysis taking too long (> {timeout}s).")
        return None
    except requests.exceptions.HTTPError as http_err:
        # Handle specific HTTP errors based on status code
        status_code = http_err.response.status_code
        print(f"Error: HTTP Error {status_code} calling Hugging Face API ({model_name}).")
        try:
            # Try to get more detailed error message from HF response body
            error_details = http_err.response.json()
            print(f"  Response Body: {error_details}")
            # Check for common specific issues
            if status_code == 401: print("  Hint: Check if HF_API_TOKEN is correct and valid.")
            if status_code == 503: print("  Hint: Model might be loading (if wait_for_model=False) or unavailable.")
            if status_code == 429: print("  Hint: Rate limit exceeded. Try again later or check HF plan.")
        except ValueError: # If response body is not JSON
            print(f"  Response Text: {http_err.response.text[:500]}...")
        return None
    except requests.exceptions.RequestException as req_err:
        # Handle other connection/request errors
        print(f"Error: Network or request issue calling Hugging Face API: {req_err}")
        return None
    except Exception as e:
        # Catch any other unexpected errors during the process
        print(f"An unexpected error occurred during Hugging Face API call:")
        traceback.print_exc()
        return None

# --- Example Usage Block ---
if __name__ == "__main__":
    print("--- Testing llm_integration/summarizer.py (Hugging Face) ---")

    if not API_TOKEN:
        print("\n--- Test Skipped ---")
        print("HF_API_TOKEN is not set in your .env file.")
    else:
        print(f"Attempting summary using model: {DEFAULT_HF_API_URL.split('/')[-1]}")
        # Use the prompt builder to create a realistic test prompt
        # Assuming prompt_builder is in the same directory or path is set
        try:
            from prompt_builder import build_llm_prompt

            # Create simple dummy data for the prompt builder
            test_analysis_data = {
                "sequence_id": "TestSeq_HF", "sequence_length": 500,
                "results": {
                    "orfs": [{"orf_id": "orf_1", "start": 50, "end": 350, "strand": "+", "length_bp": 301, "protein_sequence": "MSKGEEL..."}],
                    "motifs": [{"motif_id": "TATA_BOX_LIKE", "start": 20, "end": 27, "strand": "+", "matched_sequence": "TATAATAA"}],
                    "domains": [{"orf_id": "orf_1", "source_db": "PFAM", "accession": "PF00120", "description": "Example Domain", "start_aa": 10, "end_aa": 80, "evalue": "1e-20"}]
                }
            }
            # Build prompt using the recommended format for the default model (Mistral)
            test_prompt = build_llm_prompt(test_analysis_data, format="mistral")
            print("\nGenerated Test Prompt (Mistral Format):")
            print(test_prompt)

            # Generate the summary using the function in this file
            summary_result = generate_huggingface_summary(test_prompt)

            if summary_result:
                print("\n--- Test Summary ---")
                print(summary_result)
            else:
                print("\n--- Test Failed ---")
                print("Could not generate summary. Check API token, model URL, network, and logs.")

        except ImportError:
             print("\nError: Could not import 'build_llm_prompt'. Make sure prompt_builder.py is accessible.")
        except Exception as e:
             print("\nAn error occurred during the test setup or execution:")
             traceback.print_exc()