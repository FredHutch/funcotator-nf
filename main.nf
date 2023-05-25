#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process funcotator {
    container "${params.container__funcotator}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
        path variant_vcf
        path reference_fasta
        val ref_version
        path data_sources

    output:
        path "*.funcotated.*"

    script:
    template "funcotator.sh"

}

workflow {
    if ( params.input_vcf == false ){
        error "Must provide --input_vcf"
    }
    if ( params.outdir == false ){
        error "Must provide --outdir"
    }

    Channel
        .fromPath(
            params.input_vcf.split(',').toList()
        )
        .set {
            vcf_ch
        }
}