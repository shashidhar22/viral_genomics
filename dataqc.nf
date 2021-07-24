#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process runFastqc {
  container "$params.fastqc"
  label 'low_mem'
  publishDir "$params.out_path/fastqc/${mode}", mode : "copy"
  input:
    tuple val(sample), path(fastq_paths)
    val mode
  output:
    path "${sample}_fastqc.zip", emit: fastqc_results
  script:
    rone = fastq_paths[0]
    rtwo = fastq_paths[1]
    rthree = fastq_paths[2]
    rfour = fastq_paths[3]
    if (mode == "raw") 
      """
      zcat ${rone} ${rtwo} ${rthree} ${rfour} > ${sample}.fastq
      fastqc -q ${sample}.fastq
      """
    else if (mode == "trimmed")
      """
      cat ${rone} ${rtwo} > ${sample}.fastq
      fastqc -q ${sample}.fastq
      """
}

process runMultiQC {
  container "$params.multiqc"
  label 'low_mem'
  publishDir "$params.out_path/mutli_qc", mode : "copy"
  input:
    path multi_qc_path 
  output:
    path "multiqc_report.html", emit: mutli_qc_out
  script:
    """
    multiqc .
    """
    
}

process quast {
  container "$params.quast"
  label 'low_mem'
  publishDir "$params.out_path/quast", mode : "copy"
  input:
    each genome_assembly
    path reference_genome
  output:
    path "${sample}", emit: quast
  script:
    memory = "$task.memory" =~ /\d+/
    sample = genome_assembly.getParent().getSimpleName()
    """
    quast.py ${genome_assembly} -r ${reference_genome}  -o ${sample}
    """
    
}


/*workflow {
  fastq_paths = Channel.fromFilePairs(params.fastq_paths, size: 4)
  out_path = Channel.from(params.out_path)

  main:
    runFastqc(fastq_paths)
    runMultiQC(runFastqc.out.raw_fastqc_out.collect{"$it"})
}*/