# app.py
import streamlit as st
import os
import tempfile
import sys
import pandas as pd # For tables
import json         # For download button
import traceback    # For showing errors nicely
from dna_features_viewer import GraphicFeature, GraphicRecord # For plotting
import matplotlib.pyplot as plt # For plotting
import numpy as np              # For plotting helpers

# --- Add src directory to Python path for imports ---
app_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(app_dir, 'src')
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

# --- Import the main analysis pipeline function ---
try:
    # This function now handles all analysis including the chosen LLM summarizer
    from analysis_pipeline import run_full_analysis
except ImportError as e:
    st.error(f"Fatal Error: Could not import the analysis pipeline from 'src'.\n{e}\n"
             f"Please ensure 'src/analysis_pipeline.py' exists and 'src' is importable.\n"
             f"Current sys.path: {sys.path}")
    st.stop() # Stop app execution if core pipeline is missing


# --- Visualization Function ---
# (Include the full create_feature_plot function definition here - same as before)
def create_feature_plot(analysis_results: dict):
    """
    Generates a plot visualizing sequence features using dna_features_viewer.

    Args:
        analysis_results (dict): The dictionary containing results from run_full_analysis.

    Returns:
        matplotlib.figure.Figure | None: The generated Matplotlib figure object, or None if no data/plot possible.
    """
    if not analysis_results or 'results' not in analysis_results: return None
    seq_len = analysis_results.get("sequence_length", 0)
    if not isinstance(seq_len, int) or seq_len <= 0:
        st.warning(f"Invalid sequence length ({seq_len}), cannot generate plot.")
        return None

    all_features = []
    results = analysis_results.get("results", {})
    orfs = results.get("orfs", [])
    motifs = results.get("motifs", [])
    domains = results.get("domains", [])

    # ORF Features
    orf_color = "#4682B4" # Steel Blue
    for orf in orfs:
        try:
            start, end = int(orf.get("start")), int(orf.get("end"))
            strand_char = orf.get("strand")
            strand_val = +1 if strand_char == '+' else -1 if strand_char == '-' else 0
            prot_len = len(orf.get("protein_sequence", ""))
            label = f"ORF ({prot_len} aa)"
            all_features.append(GraphicFeature(start=start, end=end, strand=strand_val, color=orf_color, label=label))
        except (ValueError, TypeError, AttributeError): continue # Skip malformed ORFs

    # Motif Features
    motif_color = "#3CB371" # Medium Sea Green
    unique_motif_ids = sorted(list(set(m.get("motif_id", "Unknown") for m in motifs)))
    cmap = plt.cm.get_cmap('Accent', max(1, len(unique_motif_ids)))
    color_map = {label: cmap(i) for i, label in enumerate(unique_motif_ids)}
    for motif in motifs:
        try:
            start, end = int(motif.get("start")), int(motif.get("end"))
            strand_char = motif.get("strand")
            strand_val = +1 if strand_char == '+' else -1 if strand_char == '-' else 0
            label = motif.get("motif_id", "Motif")
            feature_color = color_map.get(label, motif_color)
            all_features.append(GraphicFeature(start=start, end=end, strand=strand_val, color=feature_color, label=label))
        except (ValueError, TypeError, AttributeError): continue

    # Domain Features (Mapped to DNA)
    domain_color = "#FFA07A" # Light Salmon
    orfs_dict = {orf.get("orf_id"): orf for orf in orfs if orf.get("orf_id")} # Index ORFs by ID
    unique_domain_descs = sorted(list(set(d.get("description", "Unknown") for d in domains)))
    cmap_dom = plt.cm.get_cmap('Set2', max(1, len(unique_domain_descs)))
    domain_color_map = {label: cmap_dom(i) for i, label in enumerate(unique_domain_descs)}
    for domain in domains:
        try:
            orf_id = domain.get("orf_id")
            parent_orf = orfs_dict.get(orf_id)
            if not parent_orf: continue
            start_aa, end_aa = int(domain.get("start_aa")), int(domain.get("end_aa"))
            label = domain.get("description", "Domain")
            orf_start_dna, orf_end_dna = int(parent_orf.get("start")), int(parent_orf.get("end"))
            orf_strand_char = parent_orf.get("strand")

            dna_start_offset = (start_aa - 1) * 3
            dna_end_offset = (end_aa * 3) - 1 # Inclusive end relative to start of ORF sequence

            if orf_strand_char == '+':
                domain_start_dna = orf_start_dna + dna_start_offset
                domain_end_dna = orf_start_dna + dna_end_offset
                domain_strand_val = +1
            elif orf_strand_char == '-':
                # Map relative to ORF end on forward strand
                domain_start_dna = orf_end_dna - dna_end_offset
                domain_end_dna = orf_end_dna - dna_start_offset
                domain_strand_val = -1
            else: continue

            feature_color = domain_color_map.get(label, domain_color)
            all_features.append(GraphicFeature(start=domain_start_dna, end=domain_end_dna, strand=domain_strand_val,
                                               color=feature_color, label=label, thickness=10, fontdict={'size':8}))
        except (ValueError, TypeError, AttributeError, KeyError): continue # Skip malformed domains/ORFs

    if not all_features:
        st.info("No features (ORFs, Motifs, mappable Domains) found to plot.")
        return None

    # Create plot
    try:
        all_features.sort(key=lambda f: f.start)
        record = GraphicRecord(sequence_length=seq_len, features=all_features)
        # Adjust figure width dynamically, set a minimum and maximum
        fig_width = max(8, min(25, seq_len / 100))
        fig, ax = plt.subplots(1, figsize=(fig_width, 4 + len(all_features) // 10)) # Adjust height slightly based on feature count
        record.plot(ax=ax, with_sequence_level=False) # Let plot decide width based on figsize
        ax.set_title(f"Feature Map for {analysis_results.get('sequence_id', '')}", fontsize=10)
        ax.set_xlabel("Sequence Position (bp)")
        # Try to prevent labels overlapping - might need more sophisticated logic for many features
        plt.subplots_adjust(bottom=0.2, top=0.9) # Adjust margins
        # fig.tight_layout(pad=1.5) # Can sometimes cause issues with streamlit

        return fig
    except Exception as plot_err:
         print(f"Error during plot generation: {plot_err}")
         traceback.print_exc()
         st.error(f"Could not generate plot: {plot_err}")
         return None


# --- Streamlit App Layout & Logic ---

st.set_page_config(page_title="GenoView", layout="wide", initial_sidebar_state="collapsed")

st.title("🧬 GenoView: Intelligent Sequence Annotator")
st.markdown("Upload a DNA/RNA sequence in **FASTA format** to identify Open Reading Frames (ORFs), "
            "sequence motifs, protein domains (via InterProScan), and generate an AI-powered "
            "biological summary.")

# --- Session State Initialization ---
# Use session state to store results and app state across reruns
if 'analysis_results' not in st.session_state:
    st.session_state['analysis_results'] = None
if 'analysis_running' not in st.session_state:
    st.session_state['analysis_running'] = False
if 'uploaded_filename' not in st.session_state:
     st.session_state['uploaded_filename'] = None


# --- File Uploader ---
# Use columns for better layout
col1, col2 = st.columns([3, 1]) # Give uploader more space

with col1:
    uploaded_file = st.file_uploader(
        "Choose a FASTA file (.fasta, .fa, .fna)",
        type=['fasta', 'fa', 'fna'],
        key="fasta_uploader", # Add key for stability
        # Reset state if file is removed
        on_change=lambda: st.session_state.update(analysis_results=None, analysis_running=False)
    )

with col2:
     # Add the button next to the uploader
     analyze_button = st.button(
         "Analyze Sequence",
         key="analyze_button",
         type="primary",
         disabled=st.session_state['analysis_running'] or uploaded_file is None,
         use_container_width=True # Make button fill column
     )
     st.write("") # Add some vertical space maybe?


# --- Analysis Execution Logic ---
if analyze_button and uploaded_file is not None:
    # Set state immediately
    st.session_state['analysis_running'] = True
    st.session_state['analysis_results'] = None # Clear previous results
    st.session_state['uploaded_filename'] = uploaded_file.name # Store filename


    # Save uploaded file temporarily
    # Using 'with' ensures the file handle is closed properly
    # 'delete=False' is needed on Windows to allow reading the path; manual cleanup required.
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode='wb') as tmp_file:
        try:
            tmp_file.write(uploaded_file.getvalue())
            fasta_path = tmp_file.name # Get the path to the temporary file
        except Exception as e:
             st.error(f"Error saving uploaded file temporarily: {e}")
             st.session_state['analysis_running'] = False
             st.stop() # Stop if file cannot be saved


    # Display progress using spinner
    progress_message = ('Running analysis: ORF finding, motif detection, '
                        'InterProScan domain search (can take minutes), '
                        'and LLM summarization...')
    with st.spinner(progress_message):
        analysis_start_time = time.time()
        try:
            # Run the full analysis pipeline function
            # Assumes run_full_analysis uses skip_llm=False by default now
            results = run_full_analysis(fasta_path) # Pass the path

            # Store results in session state
            st.session_state['analysis_results'] = results
            analysis_end_time = time.time()
            st.success(f'Analysis Complete! (Took {analysis_end_time - analysis_start_time:.2f} seconds)')

        except Exception as e:
            st.error("An error occurred during the analysis pipeline:")
            st.exception(e) # Display the full exception traceback in the app
            st.session_state['analysis_results'] = None # Clear results on error

        finally:
             # --- IMPORTANT: Clean up the temporary file ---
             if 'fasta_path' in locals() and os.path.exists(fasta_path):
                  try:
                      os.remove(fasta_path)
                      # print(f"Debug: Temporary file {fasta_path} removed.") # Debug log
                  except Exception as e:
                       st.warning(f"Could not remove temporary file {fasta_path}: {e}")
             # Reset running state regardless of success/failure
             st.session_state['analysis_running'] = False
             # Rerun to update the display based on new session state
             st.rerun()


# --- Display Area ---
st.divider() # Separator before results area

# Display results if available in session state
if st.session_state['analysis_results']:
    results_data = st.session_state['analysis_results']
    # Display filename associated with these results
    if st.session_state['uploaded_filename']:
         st.caption(f"Showing results for: `{st.session_state['uploaded_filename']}`")

    st.subheader("📊 Analysis Results Overview")

    # --- LLM Summary ---
    st.markdown("#### Biological Summary (AI Generated)")
    summary = results_data.get('llm_summary', 'Summary not available.')
    if summary is None:
         summary = "Summary generation was skipped or failed."
         st.warning(summary)
    elif isinstance(summary, str) and ("Error:" in summary or "Failed" in summary or "Skipped" in summary or "Could not build prompt" in summary):
        st.warning(summary)
    else:
        st.success(summary)

    # --- Display Sections using Tabs ---
    tab_plot, tab_orf, tab_motif, tab_domain, tab_raw = st.tabs([
        "📊 Feature Plot", "ORFs", "Motifs", "Domains", "Raw JSON"
    ])

    # --- Feature Plot Tab ---
    with tab_plot:
        st.markdown("##### Sequence Feature Map")
        # Add a spinner specific to plot generation as it might take a moment
        with st.spinner("Generating feature plot..."):
            try:
                feature_plot_fig = create_feature_plot(results_data)
                if feature_plot_fig:
                     # Display using st.pyplot, clear figure after use recommended
                     st.pyplot(feature_plot_fig, clear_figure=True)
                # else: message is handled inside create_feature_plot
            except Exception as plot_error:
                 st.error("An error occurred while generating the feature plot:")
                 st.exception(plot_error) # Show error details

    # --- ORF Tab ---
    with tab_orf:
        st.markdown("##### Open Reading Frames (ORFs)")
        orfs = results_data.get("results", {}).get("orfs", [])
        if orfs:
            try:
                orf_df = pd.DataFrame(orfs)
                orf_df["protein_length_aa"] = orf_df["protein_sequence"].apply(lambda x: len(x) if isinstance(x, str) else 0)
                display_columns = ["orf_id", "start", "end", "strand", "length_bp", "protein_length_aa"]
                existing_cols = [col for col in display_columns if col in orf_df.columns]
                st.dataframe(orf_df[existing_cols], use_container_width=True, hide_index=True)
            except ImportError: st.error("Pandas library needed for table display.")
            except Exception as e: st.error(f"Error displaying ORF table: {e}")
        else: st.info("No ORFs meeting the criteria were found.")

    # --- Motif Tab ---
    with tab_motif:
        st.markdown("##### Sequence Motifs")
        motifs = results_data.get("results", {}).get("motifs", [])
        if motifs:
            try:
                motif_df = pd.DataFrame(motifs)
                display_columns = ["motif_id", "start", "end", "strand", "matched_sequence"]
                existing_cols = [col for col in display_columns if col in motif_df.columns]
                st.dataframe(motif_df[existing_cols], use_container_width=True, hide_index=True)
            except ImportError: st.error("Pandas library needed for table display.")
            except Exception as e: st.error(f"Error displaying Motif table: {e}")
        else: st.info("No predefined motifs were found.")

    # --- Domain Tab ---
    with tab_domain:
        st.markdown("##### Protein Domains (InterProScan)")
        domains = results_data.get("results", {}).get("domains", [])
        if domains:
            try:
                domain_df = pd.DataFrame(domains)
                columns_to_show = ["orf_id", "source_db", "accession", "description", "start_aa", "end_aa", "evalue", "interpro_id"]
                existing_cols = [col for col in columns_to_show if col in domain_df.columns]
                st.dataframe(domain_df[existing_cols], use_container_width=True, hide_index=True)
            except ImportError: st.error("Pandas library needed for table display.")
            except Exception as e: st.error(f"Error displaying Domain table: {e}")
        else: st.info("No protein domains were found (or no ORFs suitable for search).")

    # --- Raw JSON Tab ---
    with tab_raw:
        st.markdown("##### Raw JSON Output")
        try:
            json_string = json.dumps(results_data, indent=2)
            download_filename = f"{results_data.get('sequence_id', 'genoview_results')}_analysis.json"
            st.download_button(
                label="Download Full Results (JSON)",
                data=json_string,
                file_name=download_filename,
                mime="application/json",
            )
        except ImportError: st.error("JSON library not found.")
        except Exception as e: st.error(f"Could not prepare JSON for download: {e}")

        st.json(results_data, expanded=False) # Interactive JSON viewer

# Display initial message if no file is uploaded yet
elif uploaded_file is None and not st.session_state['analysis_running']:
     st.info("Please upload a FASTA file to begin.")

# Message while analysis is running (covered by spinner mostly)
# elif st.session_state['analysis_running']:
#     st.info("Analysis is in progress...")