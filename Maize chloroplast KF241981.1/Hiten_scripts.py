import os
import subprocess
import sys
from pathlib import Path

def run_script_in_directory(directory, script_name):
    """
    Run a single script in the specified directory
    """
    # Get absolute path of current directory
    original_dir = os.getcwd()
    dir_path = Path(directory)
    
    # Check if directory exists
    if not dir_path.exists():
        print(f"Error: Directory '{directory}' does not exist")
        return False
    
    print(f"\n{'=' * 50}")
    print(f"Entering directory: {dir_path.resolve()}")
    
    try:
        # Switch to target directory
        os.chdir(dir_path)
        
        script_path = Path(script_name)
        
        # Check if script exists
        if not script_path.exists():
            print(f"Error: Script '{script_name}' does not exist")
            return False
        
        print(f"\nExecuting script: {script_path.name}")
        print(f"{'-' * 40}")
        
        # Execute script
        result = subprocess.run(
            [sys.executable, script_name],  # Use current Python interpreter
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding='utf-8'
        )
        
        # Print output
        if result.stdout:
            print(result.stdout)
        
        # Check execution result
        if result.returncode != 0:
            print(f"\n❌ Script execution failed, exit code: {result.returncode}")
            return False
        
        print(f"\n✅ Script executed successfully: {script_path.name}")
        return True
        
    finally:
        # Ensure return to original directory
        os.chdir(original_dir)
        print(f"\nReturned to parent directory: {original_dir}")

def main():
    print("=" * 50)
    print("Starting automated task execution process")
    print("=" * 50)
    
    
    print("\n>>> Step 1: Execute Circle_junction_identification script in current directory")
    if not run_script_in_directory(".", "Circle_junction_identification.PY"):
        print(f"\n❌ Circle_junction_identification execution failed, stopping execution")
        sys.exit(1)
    
    # Step 2: Execute Read_junction_read_details_file.py script in current directory
    print("\n>>> Step 2: Execute Read_junction_read_details_file.py script in current directory")
    if not run_script_in_directory(".", "Read_junction_read_details_file.py"):
        print(f"\n❌ Read_junction_read_details_file.py execution failed, stopping execution")
        sys.exit(1)
    
    # Step 3: Execute file movement in current directory 
    print("\n>>> Step 3: Execute file movement task in current directory")
    if not run_script_in_directory(".", "File_creation.py"):
        print(f"\n❌ File_creation.py execution failed, stopping execution")
        sys.exit(1)
    
    # Step 4: Execute tasks in Reference_genome directory 
    print("\n>>> Step 4: Execute tasks in Reference_genome directory")
    seq_tasks = [
        "High_confidence_junction_read_screen.py"
    ]
    for script in seq_tasks:
        if not run_script_in_directory("Reference_genome", script):
            print(f"\n❌ Sequence directory task failed at script: {script}")
            sys.exit(1)
    
    # Step 5: Execute tasks in Intron_spliced_sequences directory 
    print("\n>>> Step 5: Execute tasks in Intron_spliced_sequences directory")
    recon_tasks = [
        "High_confidence_junction_read_screen.py",  # Note: Ensure actual filename matches
        "Junction_position_conversion_from_intron-spliced_to_reference_genome-maize-chloroplast-KF241981.1.py",
        "FL_transcript_screen_from_intron-spliced_junction_reads-maize-chloroplast-KF241981.1.py"
    ]
    for script in recon_tasks:
        if not run_script_in_directory("Intron_spliced_sequences", script):
            print(f"\n❌ Reconstruct directory task failed at script: {script}")
            sys.exit(1)
    
    # Step 6: Return to parent directory and execute full-length addition 
    print("\n>>> Step 6: Execute full-length addition task in current directory")
    if not run_script_in_directory(".", "Addition_of_FL_to_reference_genome_junction_reads.py"):
        print(f"\n❌ Addition_of_FL_to_reference_genome_junction_reads.py execution failed, stopping execution")
        sys.exit(1)
    
    # Step 7: Execute tasks in Hiten_results directory
    print("\n>>> Step 7: Execute tasks in Hiten_results directory")
    merge_tasks = [
        "Position_conversion_for_transcripts_derived_from_genome_junction-maize-chloroplast-KF241981.1.py",
        "End_clustering.py",        
        "End_annotation.py"
    ]
    for script in merge_tasks:
        if not run_script_in_directory("Hiten_results", script):
            print(f"\n❌ Hiten_results directory task failed at script: {script}")
            sys.exit(1)
    
    print("\n" + "=" * 50)
    print("🎉 All tasks completed successfully!")
    print("=" * 50)

if __name__ == "__main__":
    main()