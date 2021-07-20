#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process runFastqc {
  module "Singularity"
  container "$params.container.fastqc"
  label 'low_mem'
  publishDir "$params.out_path/fastqc/${mode}", mode : "copy"
  input:
    tuple val(sample), path(fastq_paths)
    val mode
  output:
    path "*.zip", emit: fastqc_results
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
  module "Singularity"
  container "$params.container.multiqc"
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
  module "Singularity"
  container "$params.container.quast"
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