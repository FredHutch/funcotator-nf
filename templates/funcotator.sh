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
     --transcript-selection-mode "${params.transcript_selection_mode}" \
     2>&1 >> "\${OUTPUT}.log" || echo Execution Failed >> "\${OUTPUT}.log"

if [ -s "\${OUTPUT}" ]; then
    echo Compressing output
    bgzip \${OUTPUT}
fi
