import pandas as pd
import os
import concurrent.futures
import threading
import numpy as np

# Set default number of threads
DEFAULT_THREADS = 8

# Define global thread lock for synchronizing file write operations between threads
lock = threading.Lock()

def process_file(file_name, output_prefix):
    """
    Process a single file, read Excel data, split columns, and output processing results.
    """
    try:
        # Read CSV file
        df = pd.read_csv(file_name)

        # Split 'circrna_id' column using pipe '|' as separator
        df_split = df['circrna_id'].str.split('|', expand=True)
        
        # Check if split result meets expectations (5 columns)
        if df_split.shape[1] == 5:
            df_split.columns = ['circ-ID', '5‘end', '3’end', 'op/gap', 'strand']
            df_split['op/gap'] = df_split['op/gap'].astype(int)
        else:
            print(f"Error: 'circrna_id' column not properly split, column count mismatch, actual columns: {df_split.shape[1]}")
            return

        # Merge original data with split data
        df = pd.concat([df, df_split], axis=1)

        # Sort by 'read_id' column
        df.sort_values(by='read_id', inplace=True)

        # Output filename
        output_file = f'{output_prefix}-multiple-processing-result.csv'
        
        # Export first processing result: keep read_id that appear only once
        df_single = df.groupby('read_id').filter(lambda x: len(x) == 1)
        with lock:
            df_single.to_csv(output_file, index=False, mode='w', header=True)

        # Process read_id that appear twice for further analysis
        equal_records = []
        min_gap_records = []

        for read_id, group in df.groupby('read_id'):
            if len(group) == 2:
                gap_op1 = group.iloc[0]['op/gap']
                gap_op2 = group.iloc[1]['op/gap']
                
                if gap_op1 == gap_op2:
                    equal_records.append(group)
                elif abs(gap_op1 - gap_op2) >= 10:
                    # Keep record with smaller gap/op value
                    if abs(gap_op1) > abs(gap_op2):
                        min_gap_records.append(group.iloc[[1]])
                    else:
                        min_gap_records.append(group.iloc[[0]])

        # Output qualified records
        with lock:
            if equal_records:
                df_equal = pd.concat(equal_records).reset_index(drop=True)
                df_equal.to_csv(output_file, index=False, mode='a', header=False)
            if min_gap_records:
                df_min_gap = pd.concat(min_gap_records).reset_index(drop=True)
                df_min_gap.to_csv(output_file, index=False, mode='a', header=False)

        print(f'Processing completed, results saved to {output_file}')

        # Output gap_sequence
        create_gap_sequence(output_file, output_prefix)
        
        # Output M-P-junction_read_details (original gap_sequence-result)
        create_gap_sequence_result(output_file, output_prefix)
        
        # Subsequent output of duplication and circ_details
        create_duplication(output_file, output_prefix, df.columns)
        create_circ_details(f'{output_prefix}-M-P-junction_read_details.csv', output_prefix)
        
        # New: Create poly(A) file (must be called after M-P-junction_read_details is generated)
        create_poly_a_output(f'{output_prefix}-M-P-junction_read_details.csv', output_prefix)

    except Exception as e:
        print(f"Error occurred while processing file {file_name}: {e}")

def create_gap_sequence(result_file, output_prefix):
    """
    Filter rows containing nucleic acid sequences from processing results and output to gap-sequence file.
    """
    gap_sequence_file = f'{output_prefix}-gap_sequence.csv'
    
    try:
        df = pd.read_csv(result_file)
        
        # Filter condition: Sequence of non-encoded nts is not empty
        df_gap = df[df['Sequence of non-encoded nts'].notnull() & (df['Sequence of non-encoded nts'] != '')]
        
        # Output to gap_sequence file
        with lock:
            df_gap.to_csv(gap_sequence_file, index=False)
        
        print(f'gap-sequence processing completed, results saved to {gap_sequence_file}')
        
    except Exception as e:
        print(f"Error occurred while processing file {result_file}: {e}")

def create_gap_sequence_result(result_file, output_prefix):
    """
    Filter records based on specific conditions and output to M-P-junction_read_details file.
    """
    gap_sequence_result_file = f'{output_prefix}-M-P-junction_read_details.csv'
    
    try:
        df = pd.read_csv(result_file)
        
        # Filter condition 1: 'Sequence of non-encoded nts' is empty
        condition1 = df['Sequence of non-encoded nts'].isnull() | (df['Sequence of non-encoded nts'] == '')
        
        # Filter condition 2: 'Number of non-encoded nts' is 1, 2, or 3
        condition2 = df['Number of non-encoded nts'].isin([1, 2, 3])
        
        # Filter condition 3: 'Number of non-encoded nts' >= 4 and 'Percentage of A' >= 70
        condition3 = (df['Number of non-encoded nts'] >= 4) & (df['Percentage of A'] >= 70)
        
        # Combine conditions, filter rows that meet any condition
        filtered_df = df[condition1 | condition2 | condition3]
        
        # Output to M-P-junction_read_details file
        with lock:
            filtered_df.to_csv(gap_sequence_result_file, index=False)
        
        print(f'M-P-junction_read_details processing completed, results saved to {gap_sequence_result_file}')
        
    except Exception as e:
        print(f"Error occurred while processing file {result_file}: {e}")

def create_duplication(result_file, output_prefix, original_columns):
    """
    Process duplicate read_id occurrences and create crossover output to duplication file.

    Args:
        result_file (str): Input filename containing processing results.
        output_prefix (str): Output filename prefix.
        original_columns (list): Original column name list for maintaining output format consistency.
    """
    duplication_file = f'{output_prefix}-duplication.csv'
    
    try:
        df_result = pd.read_csv(result_file)
        df_duplication = df_result[df_result.duplicated('read_id', keep=False)].copy()
        df_duplication.sort_values('read_id', inplace=True)
        df_duplication['circrna_id1'] = df_duplication.groupby('read_id')['circrna_id'].shift(-1)
        df_duplication['circrna_id2'] = df_duplication.groupby('read_id')['circrna_id'].shift(1)
        for read_id, group in df_duplication.groupby('read_id'):
            if len(group) == 2:
                df_duplication.loc[group.index, 'circrna_id1'] = group['circrna_id'].iloc[::-1].values
                df_duplication.loc[group.index, 'circrna_id2'] = group['circrna_id'].values
        df_duplication.fillna('', inplace=True)
        df_duplication = df_duplication[list(original_columns) + ['circrna_id1', 'circrna_id2']]
        with lock:
            df_duplication.to_csv(duplication_file, index=False)
        print(f'duplication processing completed, results saved to {duplication_file}')
    except Exception as e:
        print(f"Error occurred while processing file {result_file}: {e}")

def create_circ_details(result_file, output_prefix):
    """
    Count circ-ID frequency from M-P-junction_read_details file and output to M-P-circ_details file.
    """
    circ_details_file = f'{output_prefix}-M-P-circ_details.csv'
    
    try:
        df_result = pd.read_csv(result_file)
        circ_freq = df_result['circrna_id'].value_counts().reset_index()
        circ_freq.columns = ['circrna-id', 'count']
        
        with lock:
            circ_freq.to_csv(circ_details_file, index=False)
        
        print(f'M-P-circ_details processing completed, results saved to {circ_details_file}')
        
    except Exception as e:
        print(f"Error occurred while processing file {result_file}: {e}")

def create_poly_a_output(input_file, output_prefix):
    """
    Filter rows where column 13 is greater than 0 from M-P-junction_read_details file
    Output to poly(A) file
    """
    poly_a_file = f'{output_prefix}-poly(A).csv'
    
    try:
        # Verify source file exists
        if not os.path.exists(input_file):
            print(f"Error: Source file {input_file} does not exist")
            return
            
        # Read M-P-junction_read_details file
        df = pd.read_csv(input_file)
        
        # Check if column count is sufficient
        if df.shape[1] < 13:
            print(f"Error: File {input_file} only has {df.shape[1]} columns, at least 13 columns required")
            return
            
        # Get column 13 name (for logging)
        target_col = df.columns[12]
        print(f"Processing column 13: {target_col}")
        
        # Filter rows where column 13 (index 12) is greater than 0
        poly_a_data = df[df.iloc[:, 12] > 0]
        
        # Output to poly(A) file
        with lock:
            poly_a_data.to_csv(poly_a_file, index=False)
        
        print(f'Successfully output {len(poly_a_data)} rows to {poly_a_file}')
        
    except Exception as e:
        print(f"Error occurred while processing poly(A) file: {str(e)}")

# ==================== First Script Main Function ====================

def analysis_data_main(max_threads=DEFAULT_THREADS):
    """
    First script main function: analysis-data with correct results, input format modification, 
    correct duplication, filename modification, warning fixes, and polyA processing
    """
    print("=" * 80)
    print("Starting first script: analysis-data with correct results, input format modification, correct duplication, filename modification, warning fixes, and polyA processing")
    print("=" * 80)
    
    # Get all files with 'junction_read_details' in current directory
    files = [f for f in os.listdir('.') if 'junction_read_details' in f and f.endswith('.csv')]

    if not files:
        print("No CSV files containing 'junction_read_details' found")
        return False

    print(f"Found {len(files)} files to process: {files}")

    # Use thread pool for batch processing
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = []
        for file in files:
            output_prefix = file.replace('junction_read_details', '')
            futures.append(executor.submit(process_file, file, output_prefix))

        # Wait for all threads to complete
        concurrent.futures.wait(futures)
    
    print("=" * 80)
    print("First script execution completed: analysis-data with correct results, input format modification, correct duplication, filename modification, warning fixes, and polyA processing")
    print("=" * 80)
    
    return True

# ==================== Second Script Functions ====================

# Optimize memory: specify appropriate data types when reading CSV
def read_csv_with_optimal_dtype(input_file):
    # Use safer data type handling
    return pd.read_csv(input_file, dtype=str, low_memory=False)

# Process circ_details data and generate merged file
def process_circ_details_files():
    files = [f for f in os.listdir() if 'circ_details' in f and (f.endswith('.csv') or f.endswith('.xlsx'))]
    
    if not files:
        print("No files containing 'circ_details' found")
        return False
        
    max_row_count = 0
    
    # Calculate total row count from relevant columns (columns 1, 3, 5, ...) in each file
    for file in files:
        df = pd.read_csv(file, dtype=str) if file.endswith('.csv') else pd.read_excel(file)
        for i in range(0, df.shape[1], 2):
            max_row_count += df.iloc[:, i].notna().sum()

    print(f"Total max row count determined from relevant columns: {max_row_count}")
    vlookup_df = pd.DataFrame()
    circ_id_values = []
    
    # Process all files, adjust DataFrame row count to previously calculated total row count
    for file in files:
        file_name = os.path.splitext(file)[0]
        column_name_1 = f'{file_name}-col1'
        column_name_2 = f'{file_name}-col2'
        
        df = pd.read_csv(file, dtype=str) if file.endswith('.csv') else pd.read_excel(file)
        print(f"Processing file: {file} with {df.shape[0]} rows")
        
        df_subset = df.iloc[:, :2].reindex(range(max_row_count)).fillna('')
        df_subset.columns = [column_name_1, column_name_2]
        circ_id_values.extend(df_subset[column_name_1].dropna().tolist())
        
        vlookup_df = pd.concat([vlookup_df, df_subset], axis=1, join='outer', ignore_index=False)
        print(f"Current vlookup_df row count after processing {file}: {vlookup_df.shape[0]}")
    
    circ_id_values = [value.strip() for value in circ_id_values if value.strip()]
    circ_id_values = list(set(circ_id_values))
    circ_id_column = circ_id_values + [''] * (len(vlookup_df) - len(circ_id_values))
    circ_id_column = [value.strip() for value in circ_id_column]
    
    while vlookup_df.shape[1] < 145:
        vlookup_df[f'Unnamed: {vlookup_df.shape[1] + 1}'] = np.nan
    
    # Add circ-ID column to DataFrame at column 36
    vlookup_df.insert(145, 'circ-ID', circ_id_column[:len(vlookup_df)])
    vlookup_df.to_csv('vlookup.csv', index=False)
    print("Merge completed, results saved to 'vlookup.csv'")
    return True

# Process VLOOKUP and filter data based on conditions
def process_vlookup_excel(input_file, vlookup_output_file, final_output_file):
    df = read_csv_with_optimal_dtype(input_file)
    
    # Safely get AL_column (circ-ID column)
    if 'circ-ID' not in df.columns:
        print("Warning: 'circ-ID' column not found, using column 73 as fallback")
        AL_column = df.iloc[:, 145] if df.shape[1] > 145 else None
    else:
        AL_column = df['circ-ID']
    
    if AL_column is None:
        print("Error: Unable to determine AL_column")
        return False
    
    # Dynamically determine result column range
    end_col = min(146, df.shape[1])
    result_df = df.iloc[:, :end_col].copy()
    
    # Perform VLOOKUP operation on first 36 column pairs
    for i in range(72):
        lookup_range_start = i * 2
        lookup_range_end = min(lookup_range_start + 2, df.shape[1])
        
        if lookup_range_end > df.shape[1]:
            print(f"Skipping column pair {i+1}, exceeds data range")
            continue
            
        lookup_df = df.iloc[:, lookup_range_start:lookup_range_end].copy()
        lookup_df.columns = ['Key', f'Value_{i+1}']
        lookup_df = lookup_df.drop_duplicates(subset='Key')
        
        # Create mapping dictionary, ensure keys are strings
        mapping_dict = dict(zip(
            lookup_df['Key'].astype(str), 
            lookup_df[f'Value_{i+1}'].astype(float)
        ))
        
        mapped_values = AL_column.astype(str).map(mapping_dict)
        result_df[f'Column_{i+1}'] = mapped_values.fillna(0)
    
    result_df.to_csv(vlookup_output_file, index=False)
    print(f'VLOOKUP results saved to {vlookup_output_file}')
    
    # Safely perform conditional filtering
    # Determine conditional filtering column range
    condition_start = min(146, result_df.shape[1])
    condition_end = min(218, result_df.shape[1])
    
    if condition_end - condition_start < 72:
        print(f"Warning: Conditional filtering columns insufficient (actual {condition_end-condition_start} columns)")
    
    aa_al_columns = result_df.iloc[:, condition_start:condition_end]
    
    # Ensure data types are numeric
    for col in aa_al_columns.columns:
        aa_al_columns[col] = pd.to_numeric(aa_al_columns[col], errors='coerce').fillna(0)
    
    row_condition = (aa_al_columns > 0).sum(axis=1) >= 2
    filtered_df = result_df[row_condition]
    filtered_df.to_csv(final_output_file, index=False)
    print(f'Conditional filtering results saved to {final_output_file}')
    
    # Process circ-ID column splitting
    circ_id_split_and_adjust(final_output_file)
    return True

# Modify column headers and delete columns before AL column (core fix section)
def modify_column_headers_and_delete_columns(input_file):
    df = pd.read_csv(input_file, dtype=str)
    original_columns = df.columns.tolist()
    
    # Ensure only operating on first 36 columns (18 column pairs)
    columns_to_copy = original_columns[:144:2]  # Take column names 0,2,4,...34 (first column of each pair)
    
    # Copy column names to Column_1 to Column_18 positions (starting from column 38)
    new_columns = original_columns.copy()
    for i, col_name in enumerate(columns_to_copy):
        target_idx = 146 + i
        if target_idx < len(new_columns):
            new_columns[target_idx] = col_name
    
    df.columns = new_columns
    
    # Delete first 36 columns (keep content starting from column 36, including circ-ID)
    if df.shape[1] > 145:
        df = df.iloc[:, 145:]
    else:
        print("Warning: Insufficient data columns, skipping column deletion")
    
    df.to_csv(input_file, index=False)
    print(f"Header modification completed and columns before AL deleted, results saved to {input_file}")
    return True

# Split circ-ID column and readjust columns (enhanced security)
def circ_id_split_and_adjust(input_file):
    try:
        # Safely read CSV, treat all columns as strings
        df = pd.read_csv(input_file, dtype=str)
        
        if 'circ-ID' not in df.columns:
            print("Error: 'circ-ID' column not found")
            # Try to find possible alternative columns
            for col in df.columns:
                if 'circ' in col.lower() and 'id' in col.lower():
                    df.rename(columns={col: 'circ-ID'}, inplace=True)
                    print(f"Found alternative column '{col}' as circ-ID")
                    break
            else:
                print("No possible circ-ID alternative columns found")
                return False
        
        # Safely split circ-ID column
        print("circ-ID column data type:", df['circ-ID'].dtype)
        print("First 5 values:", df['circ-ID'].head().values)
        
        # Split circ-ID column, keep only first 5 columns
        circ_id_split = df['circ-ID'].str.split('|', expand=True)
        
        # Ensure at least 5 columns, fill missing with NaN
        for i in range(5):
            if i >= circ_id_split.shape[1]:
                circ_id_split[i] = np.nan
        
        # Rename columns
        circ_id_split.columns = ['T', "5'end", "3'end", 'OP', 'strand']
        
        # Merge split columns and delete 'T' column
        df = pd.concat([df, circ_id_split], axis=1).drop(columns=['T'])
        
        # Force DataFrame to rewrite column names to prevent pandas automatic naming issues
        df.to_csv(input_file, index=False)
        print(f'circ-ID column splitting completed, results saved to {input_file}')
        return True
    
    except Exception as e:
        print(f"Error occurred while processing circ-ID column: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

# ==================== Second Script Main Function ====================

def screening_optimization_main():
    """
    Second script main function: screening result optimization with correct results, 
    extraction from first row, transcript splitting, correct headers, final version, 
    secondary header modification, 18 issues-1, string issues, and 36 columns
    """
    print("=" * 80)
    print("Starting second script: screening result optimization with correct results, extraction from first row, transcript splitting, correct headers, final version, secondary header modification, 18 issues-1, string issues, and 36 columns")
    print("=" * 80)
    
    try:
        # Step 1: Process circ_details files
        if not process_circ_details_files():
            print("Failed to process circ_details files, skipping subsequent steps")
            return False
        
        # Step 2: Process VLOOKUP
        if not process_vlookup_excel('vlookup.csv', 'vlookup_result.csv', 'vlookup-result-sift.csv'):
            print("VLOOKUP processing failed, skipping subsequent steps")
            return False
        
        # Step 3: Modify headers
        if not modify_column_headers_and_delete_columns('vlookup_result.csv'):
            print("Failed to modify vlookup_result.csv headers")
            return False
            
        if not modify_column_headers_and_delete_columns('vlookup-result-sift.csv'):
            print("Failed to modify vlookup-result-sift.csv headers")
            return False
        
        print("=" * 80)
        print("Second script execution completed: screening result optimization with correct results, extraction from first row, transcript splitting, correct headers, final version, secondary header modification, 18 issues-1, string issues, and 36 columns")
        print("=" * 80)
        return True
        
    except Exception as e:
        print(f"Error occurred during second script execution: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

# ==================== Combined Main Function ====================

def main(max_threads=DEFAULT_THREADS, run_second_script=True):
    """
    Combined main function, executing two scripts in strict order
    """
    
    # Step 1: Execute first script
    first_script_success = analysis_data_main(max_threads)
    
    if not first_script_success:
        print("First script execution failed, skipping second script")
        return
    
    # Step 2: Execute second script (if first script successful and user requests it)
    if run_second_script and first_script_success:
        second_script_success = screening_optimization_main()
        
        if second_script_success:
            print("=" * 80)
            print("Both scripts executed successfully!")
            print("=" * 80)
        else:
            print("Second script execution failed")
    else:
        print("Skipping second script execution")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Batch process junction_read_details files")
    parser.add_argument('-T', '--threads', type=int, default=DEFAULT_THREADS, help="Maximum number of threads, default is 8")
    parser.add_argument('--skip-second', action='store_true', help="Skip second script execution")
    
    args = parser.parse_args()

    if args.threads > 32:
        print("Maximum threads cannot exceed 32, set to 32")
        args.threads = 32

    main(max_threads=args.threads, run_second_script=not args.skip_second)