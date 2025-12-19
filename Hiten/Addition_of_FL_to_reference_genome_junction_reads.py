import os
import pandas as pd
from pathlib import Path
import shutil

def merge_csv_files():
    # Get current working directory
    current_dir = Path.cwd()
    print(f"Current working directory: {current_dir}")
    
    # Define source and target directories
    seq_dir = current_dir / "Reference_genome"
    recon_dir = current_dir / "Intron_spliced_sequences"
    merge_dir = current_dir / "Hiten_results"
    
    # Define source file paths
    sift_csv = seq_dir / "vlookup-result-sift.csv"
    intron_csv = recon_dir / "intron-FL.csv"
    
    # Define target file path
    merged_csv = merge_dir / "Filter-the-results+FL.csv"
    
    # Check if source directories and files exist
    if not seq_dir.exists():
        print(f"❌ Error: Directory '{seq_dir}' does not exist")
        return False
    
    if not recon_dir.exists():
        print(f"❌ Error: Directory '{recon_dir}' does not exist")
        return False
    
    if not sift_csv.exists():
        print(f"❌ Error: File '{sift_csv}' does not exist")
        return False
    
    if not intron_csv.exists():
        print(f"❌ Error: File '{intron_csv}' does not exist")
        return False
    
    # Create target directory if it doesn't exist
    merge_dir.mkdir(parents=True, exist_ok=True)
    print(f"✅ Created/confirmed target directory: {merge_dir}")
    
    try:
        print("\n" + "="*50)
        print(f"Processing file: {sift_csv.name}")
        print("="*50)
        
        # Read SIFT results CSV file
        sift_df = pd.read_csv(sift_csv)
        print(f"Read successfully! Rows: {len(sift_df)}, Columns: {len(sift_df.columns)}")
        
        # Write SIFT results to merged file
        sift_df.to_csv(merged_csv, index=False)
        print(f"✅ Created initial merged file: {merged_csv.name}")
        
        print("\n" + "="*50)
        print(f"Processing file: {intron_csv.name}")
        print("="*50)
        
        # Read intron FL results CSV file
        intron_df = pd.read_csv(intron_csv)
        print(f"Read successfully! Rows: {len(intron_df)}, Columns: {len(intron_df.columns)}")
        
        # Append intron FL results to merged file
        # Use mode='a' for append, header=False to skip writing column names
        intron_df.to_csv(merged_csv, mode='a', index=False, header=False)
        print(f"✅ Appended data to merged file")
        
        # Read and verify the final merged file
        merged_df = pd.read_csv(merged_csv)
        total_rows = len(merged_df)
        print(f"\nMerge completed! Total rows: {total_rows}")
        print(f"Final file location: {merged_csv}")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Error during processing: {str(e)}")
        return False

if __name__ == "__main__":
    print("="*50)
    print("Starting CSV file merge")
    print("="*50)
    
    success = merge_csv_files()
    
    print("\n" + "="*50)
    if success:
        print("✅ CSV file merge operation completed successfully!")
    else:
        print("❌ CSV file merge operation failed, please check error messages")
    print("="*50)