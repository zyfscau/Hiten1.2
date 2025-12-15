import pandas as pd
import os

# ==================== First Script: Clustering Input File Generation ====================

def generate_clustering_input_files():
    """
    Main function of the first script: Clustering input file generation with 3-5 format modification 
    for 1+18+36 issues and positive/negative strand handling
    """
    print("=" * 80)
    print("Starting first script: Clustering input file generation with 3-5 format modification")
    print("=" * 80)
    
    # Create CSV files
    for file in ["5'end.csv", "3'end.csv", "5'end-reselt.csv", "3'end-reselt.csv"]:
        open(file, 'w').close()

    # Read High_confidence_transcripts.csv
    vlookup_file = 'High_confidence_transcripts.csv'
    try:
        df_vlookup = pd.read_csv(vlookup_file)
        print(f"Successfully read file: {vlookup_file}")
    except Exception as e:
        print(f"Error reading file {vlookup_file}: {e}")
        return False

    def process_end(df_vlookup, col_A_idx, col_B_idx, output_file):
        """
        Process 5' end or 3' end data and generate final CSV file
        :param df_vlookup: Input DataFrame
        :param col_A_idx: Column A index (end position)
        :param col_B_idx: Column B index (strand position)
        :param output_file: Output filename
        """
        # Determine if processing 5' end or 3' end
        is_5prime = "5" in output_file
        
        # Select specified columns and rename - ensure first column is end, second column is strand
        df_end = df_vlookup.iloc[:, [col_A_idx, col_B_idx]].copy()
        df_end.columns = ['end', 'strand']  # Modify column names

        # Select columns 2 to 19 (indices 1 to 72), add to df_end
        df_end = pd.concat([df_end, df_vlookup.iloc[:, 1:73]], axis=1)

        # ====== Modify strand column ======
        if is_5prime:
            # 5' end processing: modify strand column
            df_end['strand'] = df_end['strand'].replace({
                '+/-': '+',
                '-/+': '+'
            })
        else:
            # 3' end processing: modify strand column
            df_end['strand'] = df_end['strand'].replace({
                '+/-': '-',
                '-/+': '+'
            })
        # =========================
        
        # Generate X column (combine end and modified strand)
        df_end['end-strand'] = df_end['end'].astype(str) + '|' + df_end['strand'].astype(str)

        # Copy end and modified strand to AN, AO
        df_end['end.1'] = df_end['end']
        df_end['strand.1'] = df_end['strand']  # Copy modified strand value

        # Copy X column
        df_end['X'] = df_end['end-strand']

        # Insert empty column New_Col
        df_end.insert(74, 'New_Col', '')

        # **Step 1: Save first 22 columns**
        df_end.iloc[:, :76].to_csv(output_file, index=False)

        # **Step 2: Remove duplicates from X column, but don't affect AN, AO columns**
        df_selected = df_end[['end.1', 'strand.1', 'X']].drop_duplicates(subset=['X'])

        # Remove all blank rows
        df_selected = df_selected.dropna()  # Remove NaN
        df_selected = df_selected[(df_selected != '').all(axis=1)]  # Remove empty strings

        # **Reset index to start from 0, avoid empty rows**
        df_selected = df_selected.reset_index(drop=True)

        # **Read first 22 columns data**
        df_result = pd.read_csv(output_file)

        # **Merge first 22 columns + unique AN, AO, X**
        df_result = pd.concat([df_result, df_selected], axis=1)

        # **Modify final output column headers**
        df_result.columns = list(df_result.columns[:-3]) + ['end.1', 'end.1', 'X']

        # **Output final CSV**
        df_result.to_csv(output_file, index=False)
        print(f"Generated file: {output_file}")

    # Process 5' end data
    process_end(df_vlookup, 73, 76, "5'end-reselt.csv")

    # Process 3' end data
    process_end(df_vlookup, 74, 76, "3'end-reselt.csv")

    # **Delete 5'end.csv and 3'end.csv**
    for file in ["5'end.csv", "3'end.csv"]:
        if os.path.exists(file):
            os.remove(file)
            print(f"Deleted temporary file: {file}")

    def create_clustering_file(reselt_file, clustering_file):
        """
        Generate clustering file based on reselt file, preserving original column names
        """
        # Read reselt file and get column names
        try:
            df_reselt = pd.read_csv(reselt_file)
            original_columns = df_reselt.columns.tolist()
            print(f"Successfully read reselt file: {reselt_file}")
        except Exception as e:
            print(f"Error reading reselt file {reselt_file}: {e}")
            return False

        # Define column indices to extract (consistent with original code logic)
        selected_indices = [
            75,   # Column 22 (end)
            2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,   # Columns 3-12
            36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,72, # Columns 13-19
            73,   # Column 20
            76,   # Column 23 (X)
            78, 77  # Column 25 (strand), Column 24
        ]

        # Get corresponding original column names
        selected_columns = [original_columns[i] for i in selected_indices]

        # Directly extract required columns (preserving original column names)
        df_clustering = df_reselt.iloc[:, selected_indices]

        # Rename last two columns (fix order)
        df_clustering.columns = selected_columns[:-2] + [original_columns[76], original_columns[77]]

        # Save file (column names consistent with reselt file)
        df_clustering.to_csv(clustering_file, index=False)
        print(f"Generated clustering file: {clustering_file}")
        return True

    # Generate clustering files
    clustering_files_created = []
    if create_clustering_file("5'end-reselt.csv", "5'end-clustering.csv"):
        clustering_files_created.append("5'end-clustering.csv")
    if create_clustering_file("3'end-reselt.csv", "3'end-clustering.csv"):
        clustering_files_created.append("3'end-clustering.csv")

    # Delete reselt files
    for file in ["5'end-reselt.csv", "3'end-reselt.csv"]:
        if os.path.exists(file):
            os.remove(file)
            print(f"Deleted temporary file: {file}")

    print("First script execution completed: Clustering input file generation with 3-5 format modification")
    print(f"Generated clustering files: {clustering_files_created}")
    
    return len(clustering_files_created) > 0

# ==================== Second Script: Clustering Function Modification ====================

def modify_clustering_files():
    """
    Main function of the second script: Clustering function modification for input files + 18+36 issues
    """
    print("\n" + "=" * 80)
    print("Starting second script: Clustering function modification for input files")
    print("=" * 80)
    
    # Get all CSV files containing 'end-clustering' in current directory
    file_list = [f for f in os.listdir('.') if 'end-clustering' in f and f.endswith('.csv')]
    
    if not file_list:
        print("No CSV files containing 'end-clustering' found")
        return False
    
    print(f"Found {len(file_list)} files containing 'end-clustering': {file_list}")

    processed_files = []
    # Process each file
    for file_name in file_list:
        # Read CSV file (add common encoding compatibility)
        try:
            df = pd.read_csv(file_name, encoding='utf-8')
            print(f"\nProcessing file: {file_name}")
        except UnicodeDecodeError:
            try:
                df = pd.read_csv(file_name, encoding='gbk')  # Try encoding commonly used in Chinese systems
                print(f"\nProcessing file: {file_name} (using GBK encoding)")
            except Exception as e:
                print(f"Unable to read file {file_name}: {e}")
                continue

        # Check DataFrame column count
        print(f"Found {len(df.columns)} columns, column names: {df.columns.tolist()}")

        # Preprocess key columns
        for col_idx in [0, 74]:  # Original A and U columns
            if col_idx < len(df.columns):
                df.iloc[:, col_idx] = df.iloc[:, col_idx].astype(str).str.strip()

        # Define input column range: B to S columns (indices 1 to 18)
        input_columns = df.columns[1:73] if len(df.columns) >= 73 else df.columns[1:]
        
        if len(input_columns) == 0:
            print(f"File {file_name} has insufficient columns, skipping processing")
            continue
            
        # Create output column names (starting from W)
        output_columns = [f'output_{chr(87 + j)}' for j in range(len(input_columns))]

        # Group by and calculate sum
        grouped_sum = df.groupby(df.columns[0])[input_columns].sum().reset_index()

        # Merge calculation results
        df = df.merge(
            grouped_sum,
            left_on=df.columns[74] if len(df.columns) > 74 else df.columns[0],
            right_on=df.columns[0],   # Grouped key
            how='left',
            suffixes=('', '_sum')
        )

        # Clean up duplicate columns
        df.drop(columns=[df.columns[0] + "_sum"], inplace=True)

        # Rename result columns
        df.rename(columns=dict(zip(input_columns, output_columns)), inplace=True)
        
        # Step 1: Delete column 21 (index 20)
        if len(df.columns) > 74:
            df.drop(df.columns[74], axis=1, inplace=True)
        
        # Step 2: Delete columns 1-19 (indices 0-18)
        if len(df.columns) >= 73:
            df.drop(df.columns[0:73], axis=1, inplace=True)

        # Save results
        result_file_name = file_name.replace('.csv', '_result.csv')
        df.to_csv(result_file_name, index=False)
        processed_files.append(result_file_name)
        print(f"Processing completed, results saved as: {result_file_name}")

    print(f"\nSecond script execution completed: Clustering function modification for input files")
    print(f"Generated result files: {processed_files}")
    
    return len(processed_files) > 0

# ==================== Combined Main Function ====================

def main():
    """
    Combined main function, executing two scripts in strict order
    """
    # Step 1: Execute first script
    first_script_success = generate_clustering_input_files()
    
    if not first_script_success:
        print("First script execution failed, skipping second script")
        return
    
    # Step 2: Execute second script
    second_script_success = modify_clustering_files()
    
    if first_script_success and second_script_success:
        print("\n" + "=" * 80)
        print("Both scripts executed successfully!")
        print("=" * 80)
    else:
        print("Some scripts failed to execute")

if __name__ == "__main__":
    main()