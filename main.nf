#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process index_vcf {
    container "${params.container__funcotator}"
    tag "${sample}"

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

process preindex_vcf {
    container "${params.container__funcotator}"
    tag "${sample}"

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

process filter_vcf {
    container "${params.container__python}"

    input:
        tuple val(sample), path("input/"), path("input/")
        path "fasta.dict"

    output:
        tuple val(sample), path("*.vcf.gz")

    script:
        template "filter_vcf.py"
}

process funcotator {
    container "${params.container__funcotator}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${sample}"

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

process bgzip {
    container "${params.container__funcotator}"
    tag "${sample}"

    input:
        tuple val(sample), path("input/")

    output:
        tuple val(sample), path("*.vcf.gz")

    """#!/bin/bash
set -e
for fp in input/*.gz; do
    echo "Block compressing \$fp"
    gunzip -c \$fp | bgzip -c > \${fp#input/}
done
"""

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

    // Index the input file to help with filtering
    preindex_vcf(variant_vcf)

    // Filter out any variants in regions which aren't
    // included in the reference genome
    filter_vcf(
        preindex_vcf.out,
        fasta_dict.out.collect()
    )

    // Convert the filtered VCF to bgzip
    bgzip(filter_vcf.out)

    // Generate the VCF index
    index_vcf(
        bgzip.out
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