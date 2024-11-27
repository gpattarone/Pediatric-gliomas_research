# Meta-analysis scripts
Useful scripts and pipeline (under construction) to collect and curate all available clinical, anatomopathological, genomic, and molecular information from pediatric patients with gliomas, including data on gene expression, mutations, genetic variations, and epigenetic profiles, from public domain databases and existing literature in order to perform corresponding comparisons with the data obtained from the patient cohort

# glioma-nextflow-pipeline
### Nextflow Script for Running BioPySRA.py
This script assumes you have Docker installed or a compatible Python environment with all dependencies from BioPySRA.py installed. The script runs BioPySRA.py for each provided accession number.

Notes:
Input File: accessions.txt should contain one SRA accession ID per line.
Outputs: Results for each SRA ID will be saved in results/biopy_sra.
Dependencies: Ensure Python3 and required libraries (prefetch, fastq-dump, etc.) are installed.

### Nextflow Script for VCF File Handling
This script assumes:

Input: Pre-existing .vcf files or CSV/TSV files to be converted into .vcf.
Tools like ANNOVAR, SnpEff, and bcftools are installed.

Notes:
Input Files: Place .vcf files in the input_vcf folder.
Output: Annotated VCFs will be saved in results/vcf_analysis/snpeff and results/vcf_analysis/annovar.
Database Paths:
Update params.snpEff_db with the correct path to your SnpEff database.
Update params.annovar_dir with the path to your ANNOVAR directory.
Parallel Processing: The | operator allows both SnpEff and ANNOVAR annotations to run in parallel.
