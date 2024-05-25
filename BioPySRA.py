import subprocess
import os

def download_sra_files(accession_number):
    print("Downloading SRA files...")
    subprocess.run(["prefetch", "-O", ".", accession_number])

def convert_to_fastq(accession_number):
    print("Converting SRA files to FastQ format...")
    subprocess.run(["fastq-dump", accession_number])

def analyze_fastq(accession_number):
    print("Analyzing FastQ files with FastQC...")
    fastq_file = accession_number + ".fastq"
    # Check if the FASTQ file exists
    if os.path.exists(fastq_file):
        subprocess.run(["fastqc", "-o", ".", fastq_file])
        # Once FastQC is finished, generate HTML report
        print("Generating HTML report...")
        subprocess.run(["fastqc", "--extract", "-f", "fastq", fastq_file])
        # Move HTML report to a separate directory
        os.makedirs("reports", exist_ok=True)
        subprocess.run(["mv", fastq_file.replace(".fastq", "_fastqc.html"), "reports"])
    else:
        print("Error: FastQ file not found.")
    print("HTML report generated in the 'reports' directory.")

def analyze_expression_rsem(accession_number):
    print("Analyzing expression with RSEM...")
    fastq_file = accession_number + ".fastq"
    if os.path.exists(fastq_file):
        # Prepare reference (assumed to be pre-built for simplicity)
        rsem_ref = "reference"  # Replace with actual RSEM reference name/path
        subprocess.run(["rsem-calculate-expression", "--paired-end", fastq_file, rsem_ref, accession_number])
    else:
        print("Error: FastQ file not found.")
    print("RSEM analysis completed.")

def extract_sequence_info(fasta_file):
    print("Extracting sequence information...")
    transcript_info_file = "transcript_information.txt"
    with open(fasta_file, "r") as fasta, open(transcript_info_file, "w") as info:
        for line in fasta:
            if line.startswith(">"):
                header = line.strip()
                sequence = next(fasta).strip()
                info.write(f"Header: {header}\n")
                info.write(f"Length: {len(sequence)}\n")
                # Add additional relevant information extraction as needed
                info.write(f"Sequence: {sequence}\n\n")
    print(f"Sequence information extracted to {transcript_info_file}")

def calculate_rpkm(accession_number):
    print("Calculating RPKM values...")
    isoform_results_file = accession_number + ".isoforms.results"
    rpkm_file = "transcript_expression.txt"
    if os.path.exists(isoform_results_file):
        with open(isoform_results_file, "r") as isoform_results, open(rpkm_file, "w") as rpkm_out:
            headers = isoform_results.readline().strip().split('\t')
            rpkm_index = headers.index("TPM")  # Replace with "RPKM" if available
            for line in isoform_results:
                columns = line.strip().split('\t')
                transcript_id = columns[0]
                rpkm_value = columns[rpkm_index]
                rpkm_out.write(f"{transcript_id}\t{rpkm_value}\n")
    else:
        print("Error: isoforms.results file not found.")
    print(f"RPKM values calculated and saved to {rpkm_file}")

if __name__ == "__main__":
    accession_number = input("Enter the SRA accession number: ")
    fasta_file = input("Enter the path to the fasta file: ")
    
    download_sra_files(accession_number)
    convert_to_fastq(accession_number)
    analyze_fastq(accession_number)
    analyze_expression_rsem(accession_number)
    extract_sequence_info(fasta_file)
    calculate_rpkm(accession_number)
