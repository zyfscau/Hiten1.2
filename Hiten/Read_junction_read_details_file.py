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
                    
                    
                    with open(xls_file, 'r', encoding='utf-8', errors='ignore') as f:
                        
                        raw_lines = [line.rstrip('\n') for line in f if line.strip()]
                    
                    if not raw_lines:
                        print(f"Warning: No valid data in {xls_file}.")
                        error_count += 1
                        continue
                    
                    
                    
                    df = pd.DataFrame([line.split('\t') for line in raw_lines])
                    
                    df.columns = df.iloc[0]
                    df = df.drop(df.index[0])
                    df = df.reset_index(drop=True)
                    
                    
                    print(f"📊 {xls_file}: {len(df)}column, {len(df.columns)}Row")
                    
                    
                    df.to_csv(
                        csv_file,
                        index=False,
                        header=True,
                        quoting=csv.QUOTE_ALL,  
                        encoding='utf-8',
                        sep=','  
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
