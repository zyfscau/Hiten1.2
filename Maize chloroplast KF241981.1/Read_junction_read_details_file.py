import os
import glob
import pandas as pd
import csv

def main():
    print("Starting conversion of junction_read_details.xls files...")
    
    patterns = ['*intron-spliced_junction_read_details', '*ref_genome_junction_read_details']
    folders = []
    for pattern in patterns:
        folders.extend(glob.glob(pattern))
    
    success_count = 0
    error_count = 0
    
    for folder in folders:
        if os.path.isdir(folder):
            xls_file = os.path.join(folder, 'junction_read_details.xls')
            csv_file = f"{folder}.csv"
            
            if os.path.isfile(xls_file):
                try:
                    # When reading, do not specify delimiter initially, let pandas auto-detect
                    # First try using sep=None and engine='python' to auto-detect delimiter
                    df = pd.read_csv(
                        xls_file,
                        sep=None,          # Auto-detect delimiter
                        header=None,       # No header (if file has header, first row will be treated as data)
                        engine='python',
                        on_bad_lines='warn',
                        dtype=str
                    )
                    
                    # Check if data is empty
                    if df.shape[0] == 0:
                        print(f"Warning: No data rows in {xls_file}.")
                    
                    # Save CSV file
                    df.to_csv(
                        csv_file,
                        index=False,
                        header=False,
                        quoting=csv.QUOTE_ALL
                    )
                    
                    print(f"✓ Successfully converted: {csv_file}")
                    success_count += 1
                    
                except Exception as e:
                    print(f"✗ Conversion failed {xls_file}: {str(e)}")
                    error_count += 1
            else:
                print(f"✗ junction_read_details.xls not found in {folder}")
                error_count += 1
    
    print(f"\nConversion completed!")
    print(f"Success: {success_count} files")
    print(f"Failed: {error_count} files")

if __name__ == "__main__":
    main()