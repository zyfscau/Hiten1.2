import pandas as pd
import glob
import os
import numpy as np

def find_files(pattern):
    """Find all matching files"""
    return glob.glob(f"*{pattern}*.csv")

# Get all cluster files and reference file
cluster_files = find_files("end-clustering_result")
refer_file = find_files("5_and_3_end_annotation")[0]  # There's only one reference file

if not cluster_files:
    raise FileNotFoundError("No clustering result files found")

# Process each cluster file
for cluster_file in cluster_files:
    print(f"\nProcessing file: {os.path.basename(cluster_file)}")
    
    # Read data
    df_cluster = pd.read_csv(cluster_file, dtype=str, keep_default_na=False)
    df_refer = pd.read_csv(refer_file, dtype=str, keep_default_na=False)
    print(f"Cluster file original row count: {len(df_cluster)}, Reference file row count: {len(df_refer)}")
    
    # Fix floating point issues
    for col in [0, 1]:
        df_cluster.iloc[:, col] = df_cluster.iloc[:, col].str.replace(r'\.0$', '', regex=True)
    
    # Generate merge key
    df_cluster.insert(74, "Merge_Key", df_cluster.iloc[:, 0] + df_cluster.iloc[:, 1])
    
    # Perform VLOOKUP to get annotation information
    refer_dict = df_refer.set_index(df_refer.columns[0])[df_refer.columns[1]].to_dict()
    df_cluster.insert(75, "Annotation", df_cluster["Merge_Key"].map(refer_dict))
    
    # Key modification: Insert complete reference columns
    # 1. Save complete reference file content as Series
    refer_col1 = df_refer.iloc[:, 0]
    refer_col2 = df_refer.iloc[:, 1]
    
    # 2. Calculate number of rows needed for expansion
    n_cluster = len(df_cluster)
    n_refer = len(df_refer)
    
    # 3. If reference file has more rows, expand cluster file
    if n_refer > n_cluster:
        # Create empty DataFrame for filling
        empty_rows = pd.DataFrame([[''] * len(df_cluster.columns)] * (n_refer - n_cluster), 
                                 columns=df_cluster.columns)
        
        # Expand cluster file
        df_cluster = pd.concat([df_cluster, empty_rows], ignore_index=True)
        print(f"Expanded cluster file to {len(df_cluster)} rows to match reference file")
    
    # 4. Insert complete reference columns
    df_cluster.insert(76, "Refer_Col1", refer_col1[:len(df_cluster)])
    df_cluster.insert(77, "Refer_Col2", refer_col2[:len(df_cluster)])
    
    # 5. If reference file has fewer rows, fill empty values
    if n_refer < len(df_cluster):
        df_cluster.loc[n_refer:, "Refer_Col1"] = ""
        df_cluster.loc[n_refer:, "Refer_Col2"] = ""
    
    # Generate output filename
    # Determine if it's 3' end or 5' end based on filename
    if "3'end-clustering_result" in cluster_file:
        output_path = "3'end_information.csv"
    elif "5'end-clustering_result" in cluster_file:
        output_path = "5'end_information.csv"
    else:
        # If neither matches, use default name
        output_path = "end_information.csv"
    
    # Save result <-- This line will now execute in all cases!
    df_cluster.to_csv(output_path, index=False, encoding='utf-8-sig')
    print(f"Generated: {output_path}, Output file row count: {len(df_cluster)}")
    print(f"Reference column 1 row count: {df_cluster['Refer_Col1'].count()}, Non-empty values: {df_cluster['Refer_Col1'].notna().sum()}")

print("\nBatch processing completed!")