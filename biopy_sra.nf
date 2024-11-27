#!/usr/bin/env nextflow

params.github_script = 'https://raw.githubusercontent.com/gpattarone/Pediatric-gliomas_research/ed98fbfcb1c6c96a2c7fa036da1e68c0fc67bcff/BioPySRA.py'
params.accessions = "accessions.txt" // Text file with one SRA accession ID per line
params.output = "./results/biopy_sra"

process run_biopy_sra {
    publishDir "${params.output}", mode: 'copy'

    input:
    path accession_id from file(params.accessions)

    output:
    path "${params.output}/${accession_id}.results"

    script:
    """
    # Download the script from GitHub
    wget -O BioPySRA.py ${params.github_script}

    # Execute the BioPySRA script with the accession ID
    python3 BioPySRA.py --accession ${accession_id} --output ${params.output}/${accession_id}.results
    """
}
