#!/usr/bin/env python3

import argparse
import os
import subprocess

def run_star(fastq_files, reference_genome, gtf, output_dir):
    star_index = os.path.join(output_dir, 'STAR_index')
    os.makedirs(star_index, exist_ok=True)
    
    # Generate STAR index
    subprocess.run([
        'STAR', '--runThreadN', '8', '--runMode', 'genomeGenerate',
        '--genomeDir', star_index, '--genomeFastaFiles', reference_genome,
        '--sjdbGTFfile', gtf, '--sjdbOverhang', '100'
    ])
    
    # Align reads
    for fastq_file in fastq_files:
        output_prefix = os.path.join(output_dir, os.path.basename(fastq_file).split('.')[0])
        subprocess.run([
            'STAR', '--runThreadN', '8', '--genomeDir', star_index,
            '--readFilesIn', fastq_file, '--outFileNamePrefix', output_prefix,
            '--outSAMtype', 'None'
        ])

def run_rsem(fastq_files, reference_genome, gtf, output_dir):
    rsem_ref = os.path.join(output_dir, 'RSEM_reference')
    os.makedirs(rsem_ref, exist_ok=True)
    
    # Prepare RSEM reference
    subprocess.run([
        'rsem-prepare-reference', '--gtf', gtf, reference_genome, rsem_ref
    ])
    
    # Run RSEM
    for fastq_file in fastq_files:
        output_prefix = os.path.join(output_dir, os.path.basename(fastq_file).split('.')[0])
        subprocess.run([
            'rsem-calculate-expression', '--single',
            '--output-genome-bam', '--no-bam-output',
            fastq_file, rsem_ref, output_prefix
        ])

def run_variant_analysis(bam_files, reference_genome, output_dir):
    for bam_file in bam_files:
        # Step 1: Call variants using bcftools
        vcf_file = bam_file.replace('.bam', '.vcf')
        subprocess.run([
            'bcftools', 'mpileup', '-f', reference_genome, bam_file,
            '|', 'bcftools', 'call', '-mv', '-Oz', '-o', f'{vcf_file}.gz'
        ], shell=True)
        
        # Step 2: Index the VCF file
        subprocess.run(['bcftools', 'index', f'{vcf_file}.gz'])
        
        # Optionally, Step 3: Convert the gzipped VCF file to plain VCF
        subprocess.run(['bcftools', 'view', f'{vcf_file}.gz', '>', vcf_file], shell=True)
def run_preprocessing(bam_files):
    for bam_file in bam_files:
        # Step 1: Sort BAM file
        sorted_bam_file = bam_file.replace('.bam', '.sorted.bam')
        subprocess.run(['samtools', 'sort', '-o', sorted_bam_file, bam_file])
        
        # Step 2: Index the sorted BAM file
        subprocess.run(['samtools', 'index', sorted_bam_file])
        
        # Step 3: Mark duplicates using Picard
        dedup_bam_file = sorted_bam_file.replace('.sorted.bam', '.dedup.bam')
        metrics_file = sorted_bam_file.replace('.sorted.bam', '.metrics.txt')
        subprocess.run([
            'picard', 'MarkDuplicates', 'INPUT=' + sorted_bam_file,
            'OUTPUT=' + dedup_bam_file, 'METRICS_FILE=' + metrics_file
        ])

def main():
    parser = argparse.ArgumentParser(description="RNA-seq analysis pipeline")
    parser.add_argument('--fastq', nargs='+', required=True, help='List of FASTQ files')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('--gtf', required=True, help='GTF annotation file')
    parser.add_argument('--output', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    fastq_files = args.fastq
    reference_genome = args.reference
    gtf = args.gtf
    output_dir = args.output
    
    os.makedirs(output_dir, exist_ok=True)
    
    run_star(fastq_files, reference_genome, gtf, output_dir)
    run_preprocessing(fastq_files)
    run_rsem(fastq_files, reference_genome, gtf, output_dir)
    
    bam_files = [os.path.join(output_dir, os.path.basename(f).split('.')[0] + '.sorted.bam') for f in fastq_files]
    run_variant_analysis(bam_files, reference_genome, output_dir)
