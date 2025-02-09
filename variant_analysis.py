#!/usr/bin/env python3

import argparse
import os
import subprocess

def run_variant_analysis(bam_files, reference_genome, output_dir):
    """
    Performs variant calling from BAM files using bcftools.
    """
    for bam_file in bam_files:
        # Prepare output filenames
        vcf_gz_file = os.path.join(output_dir, os.path.basename(bam_file).replace('.bam', '.vcf.gz'))
        vcf_file = os.path.join(output_dir, os.path.basename(bam_file).replace('.bam', '.vcf'))

        # Step 1: Call variants using bcftools
        subprocess.run(f"bcftools mpileup -f {reference_genome} {bam_file} | bcftools call -mv -Oz -o {vcf_gz_file}", shell=True)
        
        # Step 2: Index the VCF file
        subprocess.run(['bcftools', 'index', vcf_gz_file])
        
        # Step 3: Convert the gzipped VCF file to plain VCF
        subprocess.run(f"bcftools view {vcf_gz_file} > {vcf_file}", shell=True)
        
        print(f"VCF file generated: {vcf_file}")


def run_preprocessing(fastq_files, output_dir):
    """
    Preprocesses BAM files: sorting, indexing, and duplicate marking.
    """
    bam_files = []
    for fastq_file in fastq_files:
        # Define output filenames for each step
        base_name = os.path.basename(fastq_file).split('.')[0]
        bam_file = os.path.join(output_dir, f"{base_name}.bam")
        sorted_bam_file = os.path.join(output_dir, f"{base_name}.sorted.bam")
        dedup_bam_file = os.path.join(output_dir, f"{base_name}.dedup.bam")
        metrics_file = os.path.join(output_dir, f"{base_name}.metrics.txt")

        # Step 1: Sort BAM file
        subprocess.run(['samtools', 'sort', '-o', sorted_bam_file, bam_file])
        
        # Step 2: Index the sorted BAM file
        subprocess.run(['samtools', 'index', sorted_bam_file])
        
        # Step 3: Mark duplicates using Picard
        subprocess.run([
            'picard', 'MarkDuplicates',
            f'INPUT={sorted_bam_file}', f'OUTPUT={dedup_bam_file}',
            f'METRICS_FILE={metrics_file}'
        ])
        
        print(f"Preprocessed BAM file: {dedup_bam_file}")
        bam_files.append(dedup_bam_file)

    return bam_files


def main():
    """
    Main pipeline for preprocessing and variant analysis from FASTQ and BAM files.
    """
    parser = argparse.ArgumentParser(description="Pipeline for RNA-seq and variant analysis")
    parser.add_argument('--bam', nargs='+', help='List of BAM files for variant analysis')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('--output', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    reference_genome = args.reference
    output_dir = args.output
    bam_files = args.bam

    os.makedirs(output_dir, exist_ok=True)
    
    # Preprocess BAM files and perform variant analysis
    print("Starting variant analysis...")
    run_variant_analysis(bam_files, reference_genome, output_dir)
    print("Variant analysis completed.")


if __name__ == "__main__":
    main()
