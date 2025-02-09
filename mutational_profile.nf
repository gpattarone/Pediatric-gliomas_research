#!/usr/bin/env nextflow

process ParseVCF {
    input:
    path vcf_file
    output:
    path "mutational_profile.csv"

    script:
    """
    python mutational_profile.py ${vcf_file}
    """
}

workflow {
    vcf_input = file(params.vcf)
    output_csv = ParseVCF(vcf_input)
    output_csv.view()
}
