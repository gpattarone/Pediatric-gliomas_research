import vcf
import csv
from collections import Counter

def parse_vcf_to_csv(vcf_file, output_csv):
    # Open the VCF file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    variant_data = []

    # Iterate through VCF records
    for record in vcf_reader:
        gene_name = record.INFO.get('GENE', 'Unknown')  # Assuming GENE is in INFO
        variant_type = 'SNV' if len(record.REF) == 1 and len(record.ALT[0]) == 1 else 'Indel'
        if record.INFO.get('SVTYPE'):
            variant_type = record.INFO['SVTYPE']
        
        # Append gene and variant type
        variant_data.append((gene_name, variant_type))

    # Count occurrences of each (Gene, Variant_Type)
    counts = Counter(variant_data)

    # Write to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Gene', 'Variant_Type', 'Frequency']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for (gene, variant), freq in counts.items():
            writer.writerow({'Gene': gene, 'Variant_Type': variant, 'Frequency': freq})

if __name__ == "__main__":
    vcf_file = "input.vcf"  # Replace with your VCF file
    output_csv = "mutational_profile.csv"
    parse_vcf_to_csv(vcf_file, output_csv)
