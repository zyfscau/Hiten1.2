import os
import shutil
import glob

def main():
    # Create three target folders
    folders = [
        "Hiten_results",
        "Intron_spliced_sequences",
        "Reference_genome"
    ]
    
    for folder in folders:
        os.makedirs(folder, exist_ok=True)
        print(f"Created folder: {folder}")

    # Move position-corrected Excel files
    position_files = glob.glob("*intron-spliced_junction_read_details.*")
    for file in position_files:
        if file.lower().endswith(('.xlsx', '.xls', '.csv')):
            dest = os.path.join("Intron_spliced_sequences", os.path.basename(file))
            shutil.move(file, dest)
            print(f"Moved file: {file} -> Intron_spliced_sequences")

    # Move reads_details Excel files
    reads_files = glob.glob("*ref_genome_junction_read_details.*")
    for file in reads_files:
        if file.lower().endswith(('.xlsx', '.xls', '.csv')):
            dest = os.path.join("Reference_genome", os.path.basename(file))
            shutil.move(file, dest)
            print(f"Moved file: {file} -> Reference_genome")
    
    # Additional: Move 5_and_3_end_annotation files to Hiten_results
    ref_files = glob.glob("5_and_3_end_annotation.*")
    for file in ref_files:
        if file.lower().endswith(('.xlsx', '.xls', '.csv')):
            dest = os.path.join("Hiten_results", os.path.basename(file))
            shutil.move(file, dest)
            print(f"Moved file: {file} -> Hiten_results")

    # Copy all Python scripts to three folders
    py_files = glob.glob("*.py")
    for folder in folders:
        for py_file in py_files:
            # Skip the current script itself to avoid permission issues during runtime
            if os.path.abspath(py_file) == os.path.abspath(__file__):
                continue
                
            dest = os.path.join(folder, os.path.basename(py_file))
            shutil.copy2(py_file, dest)
            print(f"Copied script: {py_file} -> {folder}")

    print("\nAll operations completed!")

if __name__ == "__main__":
    main()