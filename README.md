# GenoView: Intelligent Sequence Annotator 🧬

GenoView is a web-based platform designed to streamline the annotation of DNA/RNA sequences. It integrates core bioinformatics analyses with LLM-powered interpretation to provide researchers with structural features, functional predictions, and concise biological summaries from just a FASTA file.
<div>
    <a href="https://www.loom.com/share/02d91f2f065942ddb662ebeee5484ab1">
      <p>GenoView - 14 April 2025 - Watch Video</p>
    </a>
    <a href="https://www.loom.com/share/02d91f2f065942ddb662ebeee5484ab1">
      <img style="max-width:300px;" src="https://cdn.loom.com/sessions/thumbnails/02d91f2f065942ddb662ebeee5484ab1-131f287f6c9c7f11-full-play.gif">
    </a>
  </div>
## Features ✨

*   **FASTA Upload:** Simple interface to upload single-sequence FASTA files.
*   **ORF Prediction:** Identifies Open Reading Frames on both strands using Biopython.
*   **Motif Detection:** Scans for user-defined regex-based sequence motifs (e.g., TATA-box, GC-box).
*   **Domain Identification:** Leverages the InterProScan API to find conserved protein domains in predicted ORFs.
*   **LLM Summarization:** Utilizes Large Language Models (via Hugging Face Inference API / OpenAI API - *be specific*) to generate concise biological summaries interpreting the findings.
*   **Interactive Visualization:** Plots identified features (ORFs, motifs, mapped domains) along the sequence length using `dna_features_viewer`.
*   **Tabulated Results:** Presents detailed findings for ORFs, motifs, and domains in clear tables.
*   **JSON Report:** Allows downloading the complete analysis results in a structured JSON format.

## Tech Stack 🛠️

*   **Backend:** Python 3.10
*   **Core Libraries:** Biopython, Pandas, NumPy
*   **APIs:** InterProScan REST API, Hugging Face Inference API 
*   **Web Framework:** Streamlit
*   **Plotting:** Matplotlib, dna_features_viewer
*   **Environment:** Conda

## Setup & Installation ⚙️

1.  **Clone Repository:**
    ```bash
    git clone https://github.com/your-username/genoview.git
    cd genoview
    ```
2.  **Create Conda Environment:** (Recommended: Use the provided `environment.yml`)
    ```bash
    conda env create -f environment.yml
    conda activate genoview_v2 # Or your chosen env name
    ```
    *(Alternatively, if providing requirements.txt):*
    ```bash
    # conda create --name genoview_v2 python=3.10 -c conda-forge --yes
    # conda activate genoview_v2
    # pip install -r requirements.txt
    # conda install matplotlib -c conda-forge --yes 
    ```
3.  **API Keys:** Create a `.env` file in the project root directory (`genoview/.env`). Add your API keys/email:
    ```dotenv
    # Required for Domain Finder
    EBI_EMAIL=your_email@example.com

    # Required for LLM Summary (Choose ONE based on current setup)
    HF_API_TOKEN=hf_YourHuggingFaceToken...
    # OPENAI_API_KEY=sk-YourOpenAIKey...
    ```
    **Important:** Add `.env` to your `.gitignore` file to avoid committing secrets!

## Usage 🚀

1.  **Activate Environment:**
    ```bash
    conda activate genoview_v2 # Or your env name
    ```
2.  **Run Streamlit App:**
    ```bash
    streamlit run app.py
    ```
3.  Open the local URL provided by Streamlit (usually `http://localhost:8501`) in your web browser.
4.  Upload a single-sequence FASTA file using the file uploader.
5.  Click "Analyze Sequence".
6.  View the results in the tabs and the feature plot. Download the JSON report if needed.

## Project Structure 📂
