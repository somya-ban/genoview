# app.py

import streamlit as st
import os
import tempfile # To handle temporary file saving
import sys

# --- Add src directory to path for imports ---
# This allows app.py in the root to import modules from src/
# Get the absolute path of the directory containing app.py
app_dir = os.path.dirname(os.path.abspath(__file__))
# Construct the path to the src directory
src_dir = os.path.join(app_dir, 'src')
# Add src directory to the Python path if it's not already there
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

# --- Import your analysis pipeline ---
# Now you should be able to import from src
try:
    from analysis_pipeline import run_full_analysis
except ImportError as e:
    st.error(f"Failed to import analysis pipeline: {e}. Make sure src is in the Python path.")
    # Optionally print sys.path for debugging:
    # st.write("Current sys.path:", sys.path)
    st.stop() # Stop execution if import fails


# --- Streamlit App Layout ---

st.set_page_config(page_title="GenoView", layout="wide") # Use wide layout

st.title("🧬 GenoView: Intelligent Sequence Annotator")
st.write("Upload a DNA/RNA sequence in FASTA format to analyze its features.")

# --- File Uploader ---
uploaded_file = st.file_uploader("Choose a FASTA file (.fasta, .fa, .fna)", type=['fasta', 'fa', 'fna'])

# Placeholder for results - use session state to keep results between runs
if 'analysis_results' not in st.session_state:
    st.session_state['analysis_results'] = None
if 'analysis_running' not in st.session_state:
    st.session_state['analysis_running'] = False


# --- Analysis Trigger ---
if uploaded_file is not None:
    # Display filename
    st.write(f"Uploaded file: `{uploaded_file.name}`")

    # Button to start analysis
    analyze_button = st.button("Analyze Sequence", type="primary", disabled=st.session_state['analysis_running'])

    if analyze_button:
        st.session_state['analysis_running'] = True
        st.session_state['analysis_results'] = None # Clear previous results

        # Save uploaded file temporarily to pass its path to the pipeline
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp_file:
            tmp_file.write(uploaded_file.getvalue())
            fasta_path = tmp_file.name
            st.write(f"Temporary file saved at: {fasta_path}") # Debugging line

        st.write("--- Starting Analysis Pipeline ---")
        # Use st.spinner for progress indication
        with st.spinner('Running ORF finding, motif detection, domain search (InterProScan API), and LLM summarization... This may take several minutes.'):
            try:
                # Run the full analysis function from your pipeline script
                # Make sure run_full_analysis accepts the path and returns the dict
                results = run_full_analysis(fasta_path) # Assuming skip_llm=False by default
                st.session_state['analysis_results'] = results # Store results in session state
                st.success('Analysis Complete!')

            except Exception as e:
                st.error(f"An error occurred during analysis:")
                st.exception(e) # Display the full exception traceback in the app
                st.session_state['analysis_results'] = None # Clear results on error

            finally:
                 # Clean up the temporary file
                 if os.path.exists(fasta_path):
                      os.remove(fasta_path)
                      st.write(f"Temporary file cleaned up.") # Debugging line
                 st.session_state['analysis_running'] = False # Allow button press again
                 st.rerun() # Rerun the script to update the display based on new session state


# --- Display Results (Placeholder for now) ---
if st.session_state['analysis_results']:
    results_data = st.session_state['analysis_results']
    st.divider() # Add a visual separator
    st.subheader("📊 Analysis Results")

    # --- LLM Summary ---
    st.markdown("#### Biological Summary (LLM Generated)")
    summary = results_data.get('llm_summary', 'Summary not available.')
    if "Error:" in summary or "Failed" in summary or "Skipped" in summary:
        st.warning(summary) # Display errors/skipped messages with a warning style
    else:
        st.success(summary) # Display successful summary with a success style

    # --- Display Sections using Tabs or Expanders ---
    tab_orf, tab_motif, tab_domain, tab_raw = st.tabs(["ORFs", "Motifs", "Domains", "Raw JSON"])

    # ORF Tab
    with tab_orf:
        st.markdown("##### Open Reading Frames (ORFs)")
        orfs = results_data.get("results", {}).get("orfs", [])
        if orfs:
            # Convert list of dicts to Pandas DataFrame for better display
            import pandas as pd
            # Select and rename columns for clarity
            orf_df_display = pd.DataFrame(orfs)[["start", "end", "strand", "length_bp", "protein_sequence"]]
            orf_df_display["protein_length_aa"] = orf_df_display["protein_sequence"].apply(len)
            orf_df_display = orf_df_display[["start", "end", "strand", "length_bp", "protein_length_aa"]] # Reorder
            st.dataframe(orf_df_display, use_container_width=True)
        else:
            st.info("No ORFs meeting the criteria were found.")

    # Motif Tab
    with tab_motif:
        st.markdown("##### Sequence Motifs")
        motifs = results_data.get("results", {}).get("motifs", [])
        if motifs:
            import pandas as pd
            motif_df = pd.DataFrame(motifs)
            st.dataframe(motif_df[["motif_id", "start", "end", "strand", "matched_sequence"]], use_container_width=True)
        else:
            st.info("No predefined motifs were found.")

    # Domain Tab
    with tab_domain:
        st.markdown("##### Protein Domains (InterProScan)")
        domains = results_data.get("results", {}).get("domains", [])
        if domains:
            import pandas as pd
            domain_df = pd.DataFrame(domains)
            # Select relevant columns
            columns_to_show = ["orf_id", "source_db", "accession", "description", "start_aa", "end_aa", "evalue", "interpro_id", "interpro_desc"]
            # Filter out columns that might not exist if API changes
            existing_columns = [col for col in columns_to_show if col in domain_df.columns]
            st.dataframe(domain_df[existing_columns], use_container_width=True)
        else:
            st.info("No protein domains were found (or no ORFs suitable for search).")

    # Raw JSON Tab
    with tab_raw:
        st.markdown("##### Raw JSON Output")
        # Offer a download button for the full JSON
        try:
            json_string = json.dumps(results_data, indent=2)
            # Use sequence ID or a default name for the download file
            download_filename = f"{results_data.get('sequence_id', 'genoview_results')}_analysis.json"
            st.download_button(
                label="Download Full Results (JSON)",
                data=json_string,
                file_name=download_filename,
                mime="application/json",
            )
        except Exception as e:
             st.error(f"Could not prepare JSON for download: {e}")

        st.json(results_data, expanded=False) # Still show the interactive JSON viewer


# --- Footer or other elements ---
# Keep the logic for showing "Please upload..." or "Click Analyze..."
else:
    if uploaded_file is None:
        st.info("Please upload a FASTA file to begin.")
    elif not st.session_state['analysis_running']:
         # Only show if not running and no results yet after upload
         st.write("Click 'Analyze Sequence' to start the process.")