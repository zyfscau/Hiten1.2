#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RNA-Seq data processing pipeline script (Batch mode with 4-channel parallel processing)
Supports multiple FASTQ formats: .fq/.fastq (compressed/uncompressed), _1/_2/_R1/_R2
Author: Your Name
"""

import os
import sys
import subprocess
import shutil
import datetime
from pathlib import Path
import gzip
import re
import glob
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

class RNASeqPipeline:
    def __init__(self, sample_prefix, genome_file1="NC_011033.1.fasta", 
                 genome_file2="Rice-mitochondrion-tintron-spliced-sequences-NC_011033.1.fasta", output_dir="out", processed_log="processed_samples.log"):
        self.sample_prefix = sample_prefix  
        self.input_prefix = sample_prefix    
        self.genome_file1 = genome_file1
        self.genome_file2 = genome_file2
        self.output_dir = output_dir
        self.flash_stats_file = f"{self.sample_prefix}_flash_statistics.txt"
        self.processed_log = processed_log  
        
        # Define conda environments
        self.cutadapt_env = "/home/maize/miniconda3/bin/activate cutadaptenv"
        self.rnaseq_env = "/home/maize/miniconda3/bin/activate rna-seq"
        
        # Files to clean up 
        self.files_to_remove = [
            f"{self.input_prefix}_.clean.bf.fastq",
            f"{self.input_prefix}_.clean.bf-a.fastq",
            f"{self.input_prefix}_.extendedFrags.fastq",
            f"{self.input_prefix}_.hist",
            f"{self.input_prefix}_.histogram",
            f"{self.input_prefix}_.notCombined_1.fastq",
            f"{self.input_prefix}_.notCombined_2.fastq",
            f"{self.input_prefix}_.notCombined_2.pr.fastq"
        ]

    def run_command(self, command, description=None, env=None):
        """Run a shell command with error handling"""
        if description:
            print(f"  Step: {description}")
        
        try:
            if env:
                full_command = f"source {env} && {command}"
                result = subprocess.run(full_command, shell=True, executable='/bin/bash', 
                                      capture_output=True, text=True, check=True)
            else:
                result = subprocess.run(command, shell=True, capture_output=True, 
                                      text=True, check=True)
            
            if result.stdout:
                print(f"  Output: {result.stdout.strip()}")
            return result
        except subprocess.CalledProcessError as e:
            print(f"  Error running command: {command}")
            print(f"  Error output: {e.stderr.strip()}")
            raise

    def step1_cutadapt_trimming(self):
        """Step 1: Adapter trimming using cutadapt"""
        print("  Step 1: Trimming adapters using cutadapt...")
        
        # Get the actual file paths of R1/R2 for the current sample (compatible with all formats)
        r1_input, r2_input = self._get_sample_read_files()
        if not r1_input or not r2_input:
            raise FileNotFoundError(f"Sample {self.sample_prefix} missing R1/R2 files")
        
        print(f"  Using input files: R1={r1_input}, R2={r2_input}")
        
        # Use underscores uniformly, output file names in the format: {sample_prefix}_1.clean.fastq.gz and {sample_prefix}_2.clean.fastq.gz
        command = (
            f"cutadapt -j 16 -a AGATCGGAAGAG -A AGATCGGAAGAG -m 50 --pair-filter=both "
            f"-o '{self.sample_prefix}_1.clean.fastq.gz' -p '{self.sample_prefix}_2.clean.fastq.gz' "
            f"'{r1_input}' '{r2_input}'"
        )
        
        self.run_command(command, env=self.cutadapt_env)

    def _get_sample_read_files(self):
        """Get the actual file paths of R1/R2 for the current sample (supports all formats)"""
        # Define all possible R1 suffix patterns
        r1_patterns = [
            f"{self.sample_prefix}_R1.fastq.gz",
            f"{self.sample_prefix}_R1.fq.gz",
            f"{self.sample_prefix}_1.fastq.gz",
            f"{self.sample_prefix}_1.fq.gz",
            f"{self.sample_prefix}_R1.fastq",
            f"{self.sample_prefix}_R1.fq",
            f"{self.sample_prefix}_1.fastq",
            f"{self.sample_prefix}_1.fq"
        ]
        
        # Define all possible R2 suffix patterns (one-to-one correspondence with R1)
        r2_patterns = [
            f"{self.sample_prefix}_R2.fastq.gz",
            f"{self.sample_prefix}_R2.fq.gz",
            f"{self.sample_prefix}_2.fastq.gz",
            f"{self.sample_prefix}_2.fq.gz",
            f"{self.sample_prefix}_R2.fastq",
            f"{self.sample_prefix}_R2.fq",
            f"{self.sample_prefix}_2.fastq",
            f"{self.sample_prefix}_2.fq"
        ]
        
        # Find existing R1 files (select the first matching one by priority)
        r1_file = None
        for pattern in r1_patterns:
            if os.path.exists(pattern):
                r1_file = pattern
                break
        
        # Find existing R2 files (must match the R1 format)
        r2_file = None
        if r1_file:
            # Find the index of the matched R1 pattern, then use the same index to locate the corresponding R2 pattern
            for idx, pattern in enumerate(r1_patterns):
                if pattern == r1_file:
                    r2_candidate = r2_patterns[idx]
                    if os.path.exists(r2_candidate):
                        r2_file = r2_candidate
                    break
        
        return r1_file, r2_file

    def step2_flash_merging(self):
        """Step 2: Sequence merging using FLASH"""
        print("  Step 2: Merging sequences using FLASH...")
        
        # Create FLASH statistics file
        with open(self.flash_stats_file, 'w') as f:
            f.write("FLASH Sequence Merging Statistics\n")
            f.write("========================\n")
            f.write(f"Sample: {self.sample_prefix}\n")
            f.write(f"Execution time: {datetime.datetime.now()}\n\n")
        
        
        command = (
            f"./MeCi_1.2/bin/flash '{self.sample_prefix}_2.clean.fastq.gz' "
            f"'{self.sample_prefix}_1.clean.fastq.gz' -t 32 -M 150 "
            f"-o '{self.sample_prefix}_'"
        )
        
        result = self.run_command(command, env=self.rnaseq_env)
        flash_output = result.stdout
        
        # Parse and save statistics
        self._parse_flash_statistics(flash_output)

    def _parse_flash_statistics(self, flash_output):
        """Parse FLASH output and save statistics"""
        print(f"  FLASH Output:\n{flash_output.strip()}")
        
        with open(self.flash_stats_file, 'a') as f:
            f.write("FLASH Output:\n")
            f.write("-------------\n")
            f.write(flash_output)
            f.write("\n\nKey Statistics Summary:\n")
            f.write("------------------\n")
            
            patterns = {
                'Total pairs': r'Total pairs:\s+(\d+)',
                'Merged pairs': r'Merged pairs:\s+(\d+)',
                'Unmerged pairs': r'Unmerged pairs:\s+(\d+)',
                'Percentage of combination': r'Percentage of combination:\s+([\d.]+)%'
            }
            
            stats = {}
            for key, pattern in patterns.items():
                match = re.search(pattern, flash_output)
                if match:
                    value = match.group(1)
                    stats[key] = value
                    f.write(f"{key}: {value}\n")
            
            if 'Total pairs' in stats and 'Merged pairs' in stats:
                total = int(stats['Total pairs'])
                merged = int(stats['Merged pairs'])
                merge_ratio = (merged / total) * 100 if total > 0 else 0
                f.write(f"Merge percentage: {merge_ratio:.2f}%\n")
            
            f.write("\nOutput files:\n")
            f.write("---------\n")
            for file_pattern in [f"{self.sample_prefix}_.*.fastq"]:
                for file_path in Path('.').glob(file_pattern):
                    file_size = file_path.stat().st_size
                    size_str = self._format_file_size(file_size)
                    f.write(f"- {file_path.name} ({size_str})\n")
        
        print(f"  FLASH statistics saved to: {self.flash_stats_file}")

    def _format_file_size(self, size_bytes):
        """Format file size in human readable format"""
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size_bytes < 1024.0:
                return f"{size_bytes:.1f} {unit}"
            size_bytes /= 1024.0
        return f"{size_bytes:.1f} TB"

    def step3_sequence_processing(self):
        """Step 3: Sequence reverse complement (skip the renaming step)"""
        print("  Step 3: Sequence reverse complement...")
        
        
        command = f"seqkit seq -p -r '{self.sample_prefix}_.notCombined_2.fastq' > '{self.sample_prefix}_.notCombined_2.pr.fastq'"
        self.run_command(command, env=self.rnaseq_env)

    def step4_merge_sequence_files(self):
        """Step 4: Merge all sequence files"""
        print("  Step 4: Merging sequence files...")
        
        input_files = [
            f"{self.sample_prefix}_.extendedFrags.fastq",
            f"{self.sample_prefix}_.notCombined_1.fastq",  
            f"{self.sample_prefix}_.notCombined_2.pr.fastq"  
        ]
        
        output_file = f"{self.sample_prefix}_.clean.bf.fastq"
        with open(output_file, 'w') as outfile:
            for filename in input_files:
                if os.path.exists(filename):
                    with open(filename, 'r') as infile:
                        shutil.copyfileobj(infile, outfile)
                    print(f"  Merged: {filename}")
                else:
                    print(f"  Warning: File {filename} not found, skipping")

    def step5_quality_filtering(self):
        """Step 5: Quality filtering"""
        print("  Step 5: Quality filtering...")
        self._quality_filter_fastq()

    def _quality_filter_fastq(self):
        """Quality filter FASTQ file using Python"""
        input_file = f"{self.sample_prefix}_.clean.bf.fastq"
        output_file = f"{self.sample_prefix}_.clean.bf-a.fastq"
        
        with open(input_file, 'r') as infile, gzip.open(output_file, 'wt') as outfile:
            while True:
                header = infile.readline()
                if not header:
                    break
                seq = infile.readline()
                separator = infile.readline()
                qual = infile.readline()
                
                if not qual.startswith('@'):
                    outfile.write(header + seq + separator + qual)
        
        print(f"  Quality filtering completed: {output_file}")

    def step6_final_quality_trimming(self):
        """Step 6: Final quality trimming using cutadapt"""
        print("  Step 6: Final quality trimming...")
        
        input_file = f"{self.sample_prefix}_.clean.bf-a.fastq"
        output_file = f"{self.sample_prefix}_.clean.fastq.gz"
        command = (
            f"cutadapt -j 32 -m 25 -q 20 -o '{output_file}' "
            f"'{input_file}'"
        )
        
        self.run_command(command, env=self.cutadapt_env)

    def step7_cleanup_intermediate_files(self):
        """Step 7: Clean up intermediate files (FLASH, cutadapt, seqkit)"""
        print("  Step 7: Cleaning up intermediate files...")
        
        for file_path in self.files_to_remove:
            if os.path.exists(file_path):
                try:
                    os.remove(file_path)
                    print(f"  Deleted: {file_path}")
                except Exception as e:
                    print(f"  Warning: Could not delete {file_path}: {e}")
            else:
                print(f"  Warning: File {file_path} does not exist, skipping")

    def step8_meci_analysis_genome1(self):
        """Step 8: Run MeCi analysis with first reference genome"""
        print(f"  Step 8: Running MeCi analysis - Genome 1: {self.genome_file1}")
        
        command = (
            f"perl ./MeCi_1.2/run_MeCi.pl --in1 '{self.sample_prefix}_.clean.fastq.gz' "
            f"--genome '{self.genome_file1}' --gap 100 --overlap 3 --blast_cpu 10 "
            f"--blast_evalue 2 --chunksize 200000 --max_process 32"
        )
        
        self.run_command(command, env=self.rnaseq_env)

    def step9_cleanup_temp_genome1(self):
        """Step 9: Clean up temporary files for first genome"""
        print("  Step 9: Cleaning up temp files - Genome 1...")
        
        temp_dirs = ["01.index", "02.reads", "03.split", "04.blast", "05.parse"]
        for dir_name in temp_dirs:
            dir_path = os.path.join(self.output_dir, dir_name)
            if os.path.exists(dir_path):
                shutil.rmtree(dir_path)
                print(f"  Removed: {dir_path}")

    def step10_rename_output_genome1(self):
        """Step 10: Rename output directory for first genome"""
        genome1_short = os.path.splitext(self.genome_file1)[0]  
        final_output = f"{self.sample_prefix}.mitochondrion-{genome1_short}-100~3-ref_genome_junction_read_details"
        print(f"  Step 10: Renaming output to: {final_output}")
        
        if os.path.exists(self.output_dir):
            os.rename(self.output_dir, final_output)
            return final_output
        else:
            print(f"  Warning: Output directory {self.output_dir} does not exist")
            return None

    def step11_meci_analysis_genome2(self):
        """Step 11: Run MeCi analysis with second reference genome"""
        print(f"  Step 11: Running MeCi analysis - Genome 2: {self.genome_file2}")
        
        command = (
            f"perl ./MeCi_1.2/run_MeCi.pl --in1 '{self.sample_prefix}_.clean.fastq.gz' "
            f"--genome '{self.genome_file2}' --gap 100 --overlap 3 --blast_cpu 10 "
            f"--blast_evalue 2 --chunksize 200000 --max_process 32"
        )
        
        self.run_command(command, env=self.rnaseq_env)

    def step12_cleanup_temp_genome2(self):
        """Step 12: Clean up temporary files for second genome"""
        print("  Step 12: Cleaning up temp files - Genome 2...")
        
        temp_dirs = ["01.index", "02.reads", "03.split", "04.blast", "05.parse"]
        for dir_name in temp_dirs:
            dir_path = os.path.join(self.output_dir, dir_name)
            if os.path.exists(dir_path):
                shutil.rmtree(dir_path)
                print(f"  Removed: {dir_path}")

    def step13_rename_output_genome2(self):
        """Step 13: Rename output directory for second genome"""
        genome2_short = os.path.splitext(self.genome_file2)[0]
        final_output = f"{self.sample_prefix}.mitochondrion-{genome2_short}-100~3-intron-spliced_junction_read_details"
        print(f"  Step 13: Renaming output to: {final_output}")
        
        if os.path.exists(self.output_dir):
            os.rename(self.output_dir, final_output)
            return final_output
        else:
            print(f"  Warning: Output directory {self.output_dir} does not exist")
            return None

    def run_preprocessing_steps(self):
        """Parallel processing steps: Steps 1-7 (adapter trimming to secondary quality filtering)"""
        print(f"\n{'='*60}")
        print(f"Preprocessing sample: {self.sample_prefix}")
        print(f"{'='*60}")
        
        try:
            self.step1_cutadapt_trimming()
            self.step2_flash_merging()
            self.step3_sequence_processing()
            self.step4_merge_sequence_files()
            self.step5_quality_filtering()
            self.step6_final_quality_trimming()
            self.step7_cleanup_intermediate_files()
            
            print(f"Preprocessing completed for sample: {self.sample_prefix}")
            return True
        except Exception as e:
            print(f"Error in preprocessing sample {self.sample_prefix}: {str(e)}")
            return False

    def run_meci_analysis(self):
        """Serial processing steps: Steps 8-13 (MeCi analysis)"""
        print(f"\n{'='*60}")
        print(f"MeCi analysis for sample: {self.sample_prefix}")
        print(f"{'='*60}")
        
        try:
            # First genome analysis
            self.step8_meci_analysis_genome1()
            self.step9_cleanup_temp_genome1()
            final_output1 = self.step10_rename_output_genome1()
            
            # Second genome analysis
            self.step11_meci_analysis_genome2()
            self.step12_cleanup_temp_genome2()
            final_output2 = self.step13_rename_output_genome2()
            
            # Log the sample as processed
            self._log_processed_sample()
            
            print(f"\n MeCi analysis completed for sample {self.sample_prefix}!")
            if final_output1:
                print(f"  Result 1: {final_output1}")
            if final_output2:
                print(f"  Result 2: {final_output2}")
            
            return True
        except Exception as e:
            print(f"Error in MeCi analysis for sample {self.sample_prefix}: {str(e)}")
            return False

    def _log_processed_sample(self):
        """Log samples with completed MeCi analysis to the log file (only recorded after successful analysis)"""
        with open(self.processed_log, 'a') as f:
            f.write(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\t{self.sample_prefix}\n")

    @classmethod
    def get_processed_samples(cls, log_file="processed_samples.log"):
        """Get the list of processed samples from the log file"""
        processed = set()
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                for line in f:
                    if line.strip():
                        sample = line.strip().split('\t')[1]
                        processed.add(sample)
        return processed

    @classmethod
    def find_all_sample_pairs(cls):
        """Automatically identify all paired-end file pairs in the current directory (supports 8 format combinations: 4 compressed + 4 uncompressed)"""
        # Define all possible R1 file patterns (covering all formats)
        r1_patterns = [
            "*_R1.fastq.gz",
            "*_R1.fq.gz",
            "*_1.fastq.gz",
            "*_1.fq.gz",
            "*_R1.fastq",
            "*_R1.fq",
            "*_1.fastq",
            "*_1.fq"
        ]
        
        # Define corresponding R2 file patterns (one-to-one correspondence with R1)
        r2_patterns = [
            "*_R2.fastq.gz",
            "*_R2.fq.gz",
            "*_2.fastq.gz",
            "*_2.fq.gz",
            "*_R2.fastq",
            "*_R2.fq",
            "*_2.fastq",
            "*_2.fq"
        ]
        
        sample_pairs = []
        processed_prefixes = set()  # Avoid duplicate identification of the same sample (e.g., multiple format conflicts)
        
        # Traverse by format priority (match compressed formats first, then uncompressed; match R1/R2 first, then 1/2)
        for r1_pat, r2_pat in zip(r1_patterns, r2_patterns):
            # Find all files matching the current R1 pattern
            r1_files = glob.glob(r1_pat)
            for r1_file in r1_files:
                # Extract sample prefix (remove R1/1 and suffix)
                if "_R1." in r1_file:
                    sample_prefix = r1_file.split("_R1.")[0]
                else:  # "_1." in r1_file
                    sample_prefix = r1_file.split("_1.")[0]
                
                # Skip already identified samples and extract sample prefix
                if sample_prefix in processed_prefixes:
                    continue
                
                # Construct the corresponding R2 file path
                r2_file = r2_pat.replace("*", sample_prefix)
                if os.path.exists(r2_file):
                    # Verify if the R2 file truly corresponds (avoid ambiguity caused by prefixes containing "_R1")
                    if "_R2." in r2_pat and r2_file == r1_file.replace("_R1.", "_R2."):
                        sample_pairs.append(sample_prefix)
                        processed_prefixes.add(sample_prefix)
                        print(f"Found sample pair: {sample_prefix}")
                        print(f"  R1: {r1_file}")
                        print(f"  R2: {r2_file}")
                    elif "_2." in r2_pat and r2_file == r1_file.replace("_1.", "_2."):
                        sample_pairs.append(sample_prefix)
                        processed_prefixes.add(sample_prefix)
                        print(f"Found sample pair: {sample_prefix}")
                        print(f"  R1: {r1_file}")
                        print(f"  R2: {r2_file}")
        
        # Deduplicate and sort (ensure consistent processing order)
        unique_samples = sorted(list(set(sample_pairs)))
        return unique_samples

def process_sample_preprocessing(args):
    """Function for parallel processing: run the pre-processing steps (steps 1-7) for a sample"""
    sample_prefix, genome_file1, genome_file2, output_dir, processed_log = args
    pipeline = RNASeqPipeline(
        sample_prefix=sample_prefix,
        genome_file1=genome_file1,
        genome_file2=genome_file2,
        output_dir=output_dir,
        processed_log=processed_log
    )
    success = pipeline.run_preprocessing_steps()
    return sample_prefix, success

def main():
    """Batch process all samples (4-channel parallel processing)"""
    print("="*80)
    print("RNA-Seq Batch Processing Pipeline with 4-Channel Parallel Processing")
    print("Supports: .fq/.fastq (compressed/uncompressed), _1/_2/_R1/_R2")
    print("Parallel processing for steps 1-7 (multi-sample), sequential for steps 8-13 (per sample, sample-wise ordered)")
    print("="*80)
    
    # Configure parameters (can be modified as needed)
    genome_file1 = "NC_011033.1.fasta"
    genome_file2 = "Rice-mitochondrion-tintron-spliced-sequences-NC_011033.1.fasta"
    output_dir = "out"
    processed_log = "processed_samples.log"
    num_channels = 4  
    
    # 1. Get processed samples
    processed_samples = RNASeqPipeline.get_processed_samples(processed_log)
    print(f"\n Number of processed samples:{len(processed_samples)}")
    if processed_samples:
        print(f"Processed samples: {', '.join(processed_samples)}")
    
    # 2. Find all sample pairs to be processed (supports all formats)
    all_samples = RNASeqPipeline.find_all_sample_pairs()
    print(f"\n Number of found sample pairs：{len(all_samples)}")
    
    if not all_samples:
        print("Error: No paired-end file pairs found! Supported formats：")
        print("  _1.fq.gz/_2.fq.gz, _1.fastq.gz/_2.fastq.gz")
        print("  _R1.fq.gz/_R2.fq.gz, _R1.fastq.gz/_R2.fastq.gz")
        print("  _1.fq/_2.fq, _1.fastq/_2.fastq")
        print("  _R1.fq/_R2.fq, _R1.fastq/_R2.fastq")
        sys.exit(1)
    
    # 3. Filter out unprocessed samples
    pending_samples = [s for s in all_samples if s not in processed_samples]
    print(f"\nNumber of samples to be processed: {len(pending_samples)}") 
    if pending_samples:
        print(f"Samples to be processed：{', '.join(pending_samples)}")
    else:
        print("All samples have been processed, no need to run again！")
        sys.exit(0)
    
    # 4. Parallel preprocessing steps (Steps 1-7)
    print(f"\n{'='*80}")
    print("Phase 1: Parallel Preprocessing (Steps 1-7)")
    print(f"Using {num_channels} channels for parallel processing")
    print(f"{'='*80}")
    
    
    preprocessing_args = [
        (sample, genome_file1, genome_file2, output_dir, processed_log) 
        for sample in pending_samples
    ]
    
    # Use process pool for parallel processing
    preprocessing_success = {}
    with ProcessPoolExecutor(max_workers=num_channels) as executor:
        # Submit all tasks
        future_to_sample = {
            executor.submit(process_sample_preprocessing, args): args[0] 
            for args in preprocessing_args
        }
        
        # Collect results
        for future in as_completed(future_to_sample):
            sample = future_to_sample[future]
            try:
                sample_prefix, success = future.result()
                preprocessing_success[sample_prefix] = success
                status = "SUCCESS" if success else "FAILED"
                print(f"Preprocessing {sample_prefix}: {status}")
            except Exception as e:
                print(f"Preprocessing {sample} generated an exception: {e}")
                preprocessing_success[sample] = False
    
    # 5. Serial MeCi analysis steps (Steps 8-13)
    print(f"\n{'='*80}")
    print("Phase 2: Sequential MeCi Analysis (Steps 8-13)")
    print("Running sequentially: both genome analyses per sample are serial, and samples are processed one by one")
    print(f"{'='*80}")
    
    meci_success = {}
    for sample_prefix in pending_samples:
        if preprocessing_success.get(sample_prefix, False):
            print(f"\n Starting MeCi analysis for: {sample_prefix}")
            pipeline = RNASeqPipeline(
                sample_prefix=sample_prefix,
                genome_file1=genome_file1,
                genome_file2=genome_file2,
                output_dir=output_dir,
                processed_log=processed_log
            )
            success = pipeline.run_meci_analysis()
            meci_success[sample_prefix] = success
        else:
            print(f"\n Skipping MeCi analysis for {sample_prefix} (preprocessing failed)")
            meci_success[sample_prefix] = False
    
    # 6. Output batch processing summary
    print(f"\n{'='*80}")
    print("Batch Processing Summary")
    print(f"{'='*80}")
    
    total = len(pending_samples)
    preprocessing_success_count = sum(1 for s in preprocessing_success.values() if s)
    meci_success_count = sum(1 for s in meci_success.values() if s)
    
    print(f"Total number of samples：{total}")
    print(f"Processing successful：{preprocessing_success_count}")
    print(f"MeCi analysis successful：{meci_success_count}")
    print(f"Complete success：{meci_success_count}")
    
    # Output failed samples
    failed_preprocessing = [s for s in pending_samples if not preprocessing_success.get(s, False)]
    failed_meci = [s for s in pending_samples if preprocessing_success.get(s, False) and not meci_success.get(s, False)]
    
    if failed_preprocessing:
        print(f"\n Preprocessing failed samples：{', '.join(failed_preprocessing)}")
    if failed_meci:
        print(f"MeCi analysis failed samples：{', '.join(failed_meci)}")
    
    print(f"\nProcessed samples log：{processed_log}")
    print("="*80)

if __name__ == "__main__":
    main()