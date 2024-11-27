#!/usr/bin/env nextflow

params.vcf_dir = "./input_vcf"
params.output = "./results/vcf_analysis"
params.snpEff_db = "path_to_snpeff_database"
params.annovar_dir = "/path/to/annovar"

process load_vcf {
    publishDir "${params.output}/loaded", mode: 'copy'

    input:
    path vcf_file from file("${params.vcf_dir}/*.vcf")

    output:
    path "${params.output}/loaded/${vcf_file.getBaseName()}_loaded.vcf"

    script:
    """
    # Validate and standardize VCF file
    bcftools view ${vcf_file} -o ${params.output}/loaded/${vcf_file.getBaseName()}_loaded.vcf
    """
}

process annotate_with_snpeff {
    publishDir "${params.output}/snpeff", mode: 'copy'

    input:
    path vcf_file from load_vcf.out

    output:
    path "${params.output}/snpeff/${vcf_file.getBaseName()}_annotated.vcf"

    script:
    """
    # Annotate VCF file with SnpEff
    java -Xmx4g -jar snpEff.jar ${params.snpEff_db} ${vcf_file} > ${params.output}/snpeff/${vcf_file.getBaseName()}_annotated.vcf
    """
}

process annotate_with_annovar {
    publishDir "${params.output}/annovar", mode: 'copy'

    input:
    path vcf_file from load_vcf.out

    output:
    path "${params.output}/annovar/${vcf_file.getBaseName()}_annovar.txt"

    script:
    """
    # Annotate VCF file with ANNOVAR
    ${params.annovar_dir}/table_annovar.pl ${vcf_file} ${params.annovar_dir}/humandb/ -buildver hg19 -out ${params.output}/annovar/${vcf_file.getBaseName()}_annovar -remove -protocol refGene -operation g
    """
}

workflow {
    load_vcf | annotate_with_snpeff & annotate_with_annovar
}
