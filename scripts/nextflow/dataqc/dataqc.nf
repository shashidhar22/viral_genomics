#!/usr/bin/env nextflow


nextflow.enable.dsl=2



process runFastqc {
  module 'FastQC'
  label 'low_mem'
  publishDir "$params.out_path/raw_fastqc", mode : "copy"
  input:
    tuple val(sample), path(fastq_paths)
  output:
    path "*.zip", emit: raw_fastqc_out
  script:
    rone = fastq_paths[0]
    rtwo = fastq_paths[1]
    rthree = fastq_paths[2]
    rfour = fastq_paths[3]
    """
    zcat ${rone} ${rtwo} ${rthree} ${rfour} > ${sample}.fastq
    fastqc ${sample}.fastq
    """
}

process runMultiQC {
  module 'MultiQC'
  label 'low_mem'
  publishDir "$params.out_path/mutli_qc", mode : "copy"
  input:
    path raw_fastqc_paths 
  output:
    path "multiqc_report.html", emit: mutli_qc_out
  script:
    """
    multiqc .
    """
    
}

workflow {
  fastq_paths = Channel.fromFilePairs(params.fastq_paths, size: 4)
  out_path = Channel.from(params.out_path)

  main:
    runFastqc(fastq_paths)
    runMultiQC(runFastqc.out.raw_fastqc_out.collect{"$it"})
}