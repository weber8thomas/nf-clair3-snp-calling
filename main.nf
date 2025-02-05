#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters with default values
params.samplesheet = null         // Path to samplesheet CSV
params.output_dir = ''            // Output directory path
params.model_name = 'r941_prom_hac_g360+g422'  // Model name
params.reference = ''             // Reference genome path
params.platform = 'ont'           // Sequencing platform (ont, hifi, or ilmn)
// params.tmpdir = '/tmp'  // Temporary directory

// Print pipeline info
log.info """
====================================
Clair3 SNP Calling Pipeline
====================================
samplesheet  : ${params.samplesheet}
reference    : ${params.reference}
output dir   : ${params.output_dir}
threads      : ${params.threads}
model        : ${params.model_name}
platform     : ${params.platform}
"""

// Validate required parameters
if (!params.samplesheet) {
    error "Samplesheet parameter '--samplesheet' is required"
}
if (!params.output_dir) {
    error "Output directory parameter '--output_dir' is required"
}
if (!params.reference) {
    error "Reference genome parameter '--reference' is required"
}

// Create channel from samplesheet
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> 
        def patient_id = row.patient
        def cram = file(row.cram)
        def crai = file(row.crai)
        
        // Validate files exist
        if (!cram.exists()) {
            error "CRAM file not found: ${cram}"
        }
        if (!crai.exists()) {
            error "CRAI file not found: ${crai}"
        }
        
        return [ patient_id, cram, crai ]
    }
    .set { samples_ch }

// Process to run Clair3
process runClair3 {
    tag "${patient_id}"  // Add tag for better job tracking
    container 'hkubal/clair3:latest'
    publishDir "${params.output_dir}/${patient_id}", mode: 'copy'

    
    input:
    tuple val(patient_id), path(cram), path(crai)
    path reference
    path reference_fai

    output:
    path "*"

    script:
    """
    echo "Starting Clair3 for ${patient_id}"
    echo "CRAM: ${cram}"
    echo "Reference: ${reference}"
    echo "Using TMPDIR: \$TMPDIR"
    
    /opt/bin/run_clair3.sh \
        --bam_fn=${cram} \
        --ref_fn=${reference} \
        --threads=${params.cpus} \
        --platform="${params.platform}" \
        --model_path="/opt/models/${params.model_name}" \
        --output=.
    
    echo "Finished Clair3 for ${patient_id}"
    """
}

workflow {
    // Get reference and index files
    reference = file(params.reference)
    reference_fai = file("${params.reference}.fai")
    
    if (!reference.exists()) {
        error "Reference file not found: ${params.reference}"
    }
    if (!reference_fai.exists()) {
        error "Reference index file not found: ${params.reference}.fai"
    }
    
    // Run Clair3
    runClair3(samples_ch, reference, reference_fai)
}