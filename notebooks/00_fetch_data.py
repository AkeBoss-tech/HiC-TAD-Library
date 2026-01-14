import os
import subprocess
import sys

def fetch_data():
    # Setup paths
    # Assumes this script is in notebooks/ or similar, so we find root
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, 'data', 'raw')
    os.makedirs(data_dir, exist_ok=True)
    
    print(f"Data directory: {data_dir}")

    # 1. Download HFF Micro-C Test Data
    # URL found from cooltools OSF repository
    hff_url = "https://osf.io/3h9js/download"
    hff_file = os.path.join(data_dir, "test.mcool")
    
    print("\n[1/2] Checking for 'HFF_MicroC' test data...")
    if os.path.exists(hff_file):
        print(f"Test file already exists at: {hff_file}")
    else:
        print("Downloading test data via wget...")
        try:
            subprocess.run(["wget", hff_url, "-O", hff_file], check=True)
            print(f"Success! Test file saved to: {hff_file}")
        except Exception as e:
            print(f"Error downloading test data: {e}")

    # 2. Download GM12878 (Gold Standard) - Optional
    gm12878_url = "https://data.4dnucleome.org/files-processed/4DNFI1UEG1O1/@@download/4DNFI1UEG1O1.mcool"
    gm12878_file = os.path.join(data_dir, "4DNFI1UEG1O1.mcool")
    
    print("\n[2/2] Checking for GM12878 Gold Standard dataset...")
    if os.path.exists(gm12878_file):
        print(f"GM12878 file already exists at: {gm12878_file}")
    else:
        print(f"File not found. To download (5GB+), use the command below in your terminal:")
        print(f"wget {gm12878_url} -O {gm12878_file}")

if __name__ == "__main__":
    fetch_data()
