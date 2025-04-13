from Bio.Seq import Seq
import requests # Test requests import too
import os      # Test standard library import

print("Attempting to import core libraries...")

try:
    # Test Biopython
    my_seq = Seq("AGTACACTGGT")
    print(f"Biopython imported successfully! Sequence: {my_seq}")
    print(f"Biopython Sequence Complement: {my_seq.complement()}")

    # Test requests
    print(f"Requests imported successfully! Version: {requests.__version__}")

    # Test python-dotenv (just import for now)
    from dotenv import load_dotenv
    print("python-dotenv imported successfully!")

    print("\nSetup Test Successful! Core libraries are installed and importable.")

except ImportError as e:
    print(f"\nError importing libraries: {e}")
    print("Please check your environment activation and installation steps.")
except Exception as e:
    print(f"\nAn unexpected error occurred: {e}")