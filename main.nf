#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process funcotator {
    container "${params.container__funcotator}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(variant_vcf)
        path reference_fasta
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

    data_sources = file(
        "${params.data_sources}",
        type: 'dir',
        checkIfExists: true
    )

    funcotator(
        variant_vcf,
        reference_fasta,
        data_sources
    )
}