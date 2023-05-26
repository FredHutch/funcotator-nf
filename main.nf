#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process fasta_dict {
    container "${params.container__funcotator}"

    input:
        path reference_fasta

    output:
        path "*.dict"

"""#!/bin/bash
set -e
gatk CreateSequenceDictionary -R "${reference_fasta}"
echo Done
ls -lahtr
"""

}

process index_vcf {
    container "${params.container__funcotator}"

    input:
        tuple val(sample), path(variant_vcf)

    output:
        tuple val(sample), path(variant_vcf), path("*.tbi")

"""#!/bin/bash
set -e
gatk IndexFeatureFile -I "${variant_vcf}"
echo Done
ls -lahtr
"""

}

process funcotator {
    container "${params.container__funcotator}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(variant_vcf), path(variant_vcf_index)
        path reference_fasta
        path reference_fasta_index
        path reference_fasta_dict
        path data_sources

    output:
        path "*.funcotated.*"

    script:
    template "funcotator.sh"

}

workflow {
    if ( params.samplesheet == false ){
        error "Must provide --samplesheet"
    }
    if ( params.outdir == false ){
        error "Must provide --outdir"
    }
    if ( params.ref_version == false ){
        error "Must provide --ref_version"
    }
    if ( params.data_sources == false ){
        error "Must provide --data_sources"
    }

    // Parse the sample sheet
    Channel
        .from(
            file(
                params.samplesheet,
                checkIfExists: true
            )
        ).splitCsv(
            header: true
        ).map {
            r -> [
                r["sample"], 
                file(r["vcf"], checkIfExists: true)
            ]
        }
        .set {
            variant_vcf
        }

    reference_fasta = file(
        "${params.reference_fasta}",
        checkIfExists: true
    )

    reference_fasta_index = file(
        "${params.reference_fasta}.fai",
        checkIfExists: true
    )

    fasta_dict(reference_fasta)

    data_sources = file(
        "${params.data_sources}",
        type: 'dir',
        checkIfExists: true
    )

    // Generate the VCF index
    index_vcf(
        variant_vcf
    )

    // Run Funcotator
    funcotator(
        index_vcf.out,
        reference_fasta,
        reference_fasta_index,
        fasta_dict.out,
        data_sources
    )
}