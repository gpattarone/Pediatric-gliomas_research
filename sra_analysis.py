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

if __name__ == "__main__":
    accession_number = input("Enter the SRA accession number: ")
    download_sra_files(accession_number)
    convert_to_fastq(accession_number)
    analyze_fastq(accession_number)
