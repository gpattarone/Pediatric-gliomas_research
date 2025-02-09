# README

## Project Description
This pipeline is designed to analyze genomic and transcriptomic data from pediatric gliomas, process VCF files to generate mutational profiles, load therapy meta-analysis data into PostgreSQL, and organize genetic alterations in MongoDB.

## Files and Components

### 1. **BioPySRA**
- **Description**: Python script developed to automate the download of sequencing data from the public SRA database.
- **Key Functions**:
  - Downloads SRA data using accession numbers.
  - Converts SRA format files to FASTQ format.
  - Performs initial preprocessing to ensure data quality.
- **Input**: List of SRA accession numbers.
- **Output**: FASTQ files ready for downstream analysis.

---

### 2. **mutational_profile.py**
- **Description**: Python script that processes VCF files to extract mutational profiles and stores them in a CSV file.
- **Key Functions**:
  - Reads VCF files and extracts information about genes and variant types (SNV, Indels, CNVs).
  - Generates the `mutational_profile.csv` file with the following columns:
    - `Gene`: Name of the affected gene.
    - `Variant_Type`: Type of variant (SNV, Indel, CNV).
    - `Frequency`: Frequency of the variant.
- **Input**: VCF file.
- **Output**: `mutational_profile.csv` file.

---

### 3. **mutational_profile.nf**
- **Description**: Nextflow module to execute the `mutational_profile.py` script as part of the pipeline.
- **Key Functions**:
  - Processes a VCF file using `mutational_profile.py`.
  - Produces the `mutational_profile.csv` file.
- **Inputs**:
  - VCF file as input.
- **Outputs**:
  - `mutational_profile.csv` file.

---

### 4. **meta_analysis_therapies.csv**
- **Description**: CSV file containing a meta-analysis of targeted therapies for pediatric gliomas.
- **File Format**:
  - Columns:
    - `Therapy`: Name of the therapy.
    - `Target_Group`: Description of the target group (e.g., specific mutation or genetic alteration).
    - `Clinical_Evidence`: Summary of clinical evidence supporting the therapy.
    - `Mechanism`: Mechanism of action of the therapy.

---

### 5. **meta_analysis_therapies.nf**
- **Description**: Nextflow module to load therapy meta-analysis data into a PostgreSQL database.
- **Key Functions**:
  - Reads `meta_analysis_therapies.csv`.
  - Verifies the required column structure:
    - `Therapy`, `Target_Group`, `Clinical_Evidence`, `Mechanism`.
  - Stores the data in a PostgreSQL database.
- **Inputs**:
  - `meta_analysis_therapies.csv`.
- **Outputs**:
  - Therapy data stored in PostgreSQL.

---

### 6. **MongoDB Genetic Alterations JSON**
- **Description**: JSON file designed to store information about genetic alterations in pediatric gliomas in a MongoDB database.
- **Key Data Stored**:
  - Gene names.
  - Associated mutations.
  - Therapeutic targets.
  - Pathways and classifications based on WHO CNS5.
- **Purpose**: Enables fast and flexible querying of genetic alterations and their therapeutic implications.

---

## Usage
1. **Preprocessing**:
   - Use `BioPySRA` to download and preprocess sequencing data.
2. **Mutational Profile**:
   - Run `mutational_profile.nf` to process VCF files and generate the `mutational_profile.csv`.
3. **Therapy Meta-Analysis**:
   - Use `meta_analysis_therapies.nf` to load therapy data from `meta_analysis_therapies.csv` into PostgreSQL.
4. **Data Storage**:
   - Store genetic alteration data in MongoDB using the provided JSON structure for enhanced querying and analytics.

---

This repository provides an end-to-end pipeline for processing pediatric glioma data, ensuring reproducibility and enabling integration of genomic insights with therapeutic strategies.
