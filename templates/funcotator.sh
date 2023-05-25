#!/bin/bash

set -e

OUTPUT="${sample}.funcotated.${params.output_file_format.toLowerCase()}"
echo "Writing to \${OUTPUT}"

gatk Funcotator \
     --variant "${variant_vcf}" \
     --reference "${reference_fasta}" \
     --ref-version "${params.ref_version}" \
     --data-sources-path "${data_sources}" \
     --output "\${OUTPUT}" \
     --output-file-format "${params.output_file_format}" \
     --ignore-filtered-variants "${params.ignore_filtered_variants}" \
     --transcript-selection-mode "${params.transcript_selection_mode}" \
     --allow-hg19-gencode-b37-contig-matching ${params.allow_hg19_gencode_b37_contig_matching}

echo Compressing output
gzip \${OUTPUT}
