#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process runFastqc {
  container "$params.fastqc"
  errorStrategy 'retry'
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
      fastqc -q ${sample}.fastq > stdout.log 2> stderr.log
      """
}

process runMultiQC {
  container "$params.multiqc"
  errorStrategy 'retry'
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

process prep_quast {
  container "$params.quast"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/contigs", mode : "copy"
  input:
    tuple val(sample), path(genome_assembly)
  output:
    path("${sample}.fasta"), emit: quast_inputs
  script:
    memory = "$task.memory" =~ /\d+/
    """
    cp ${genome_assembly}/contigs.fasta ${sample}.fasta
    """
    
}


process quast {
  container "$params.quast"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/quast", mode : "copy"
  input:
    val mode
    path genome_assemblies
    path reference_genome
  output:
    path "${mode}", emit: quast
  script:
    memory = "$task.memory" =~ /\d+/
    """
    quast.py -r ${reference_genome}/genome.fa  -o ${mode} *.fasta > stdout.log 2> stderr.log
    """
    
}

process getFastq {
  module "SRA-Toolkit"
  errorStrategy 'retry'
  label 'low_mem'
  publishDir "$params.public_fastq", mode : "copy"
  input:
    each sra_number
  output:
    tuple val(sra_number), path("*.fastq") , emit: public_fastq
  script:
    """
    fasterq-dump ${sra_number} --split-3 -m 4GB -e 20 -L 0
    """
}

/*workflow {
  fastq_paths = Channel.fromFilePairs(params.fastq_paths, size: 4)
  out_path = Channel.from(params.out_path)

  main:
    runFastqc(fastq_paths)
    runMultiQC(runFastqc.out.raw_fastqc_out.collect{"$it"})
}*/