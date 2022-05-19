#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Run FASTQC on raw and trimmed FASTQ files
//
// PROCESS INPUTS:
// - Tuple with sample name and path to FASTQ files
//  (Source: fastq_paths/bbduk.out.trimmed_reads)
// - Status of FASTQ files ("trimmed", "raw")
//  (Source: "trimmed"/"raw")
// PROCESS OUTPUTS:
// - FASTQC report for FASTQ files
// EMITTED CHANNELS:
// - fastqc_results :- FASTQC reports
//
// NOTE: 
//
// TODO: 
process runFastQC {
  container "$params.fastqc"
  errorStrategy 'retry'
  label 'low_mem'
  time  '1d'
  publishDir "$params.out_path/fastqc/${mode}/results", mode : "copy", pattern : "*.zip"
  publishDir "$params.out_path/fastqc/${mode}/report", mode : "copy", pattern : "*.html"
  input:
    tuple val(sample), path(fastq_paths)
    val mode
  output:
    path "*.zip", emit: fastqc_results
    path "*.html", emit: fastqc_reports
  script:
    rone = fastq_paths[0]
    rtwo = fastq_paths[1]
    rthree = fastq_paths[2]
    rfour = fastq_paths[3]
    extension = rone.getExtension()
    """
    fastqc -t 200 --noextract -q  -o ./ *.f*q*  > stdout.txt 2> stderr.txt
    """
}

// Run MULTIQC
//
// PROCESS INPUTS:
// - Paths to run MultiQC on
//  (Source: )
// PROCESS OUTPUTS:
// - FASTQC report for FASTQ files
// EMITTED CHANNELS:
// - fastqc_results :- FASTQC reports
//
// NOTE: 
//
// TODO: 
process runMultiQC {
  container "$params.multiqc"
  errorStrategy 'retry'
  label 'low_mem'
  publishDir "$params.out_path/mutli_qc", mode : "copy"
  input:
    path raw_fastqc_out 
    path trimmed_fastqc_out
    path bbduk_stats
    path host_bowtie_stats
    path ebv_bowtie_stats
    path ebv_fastqc_out
    path prokka_spadeess
    path prokka_unicycler 
    path config
  output:
    path "multiqc_report.html", emit: mutli_qc_out
  script:
    """
    multiqc -c ${config} .
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
    path genome_assemblies
    path reference_genome
    val mode
  output:
    path "${mode}", emit: quast
  script:
    memory = "$task.memory" =~ /\d+/
    """
    quast.py -r ${reference_genome}/genome.fa --circos --glimmer -o ${mode} *.fasta > stdout.log 2> stderr.log
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

process generateFastq {
  container "$params.art"
  errorStrategy 'retry'
  label 'low_mem'
  publishDir "$params.out_path/simulated_reads", mode : "copy"
  input:
    each index
    path reference_path
  output:
    tuple val(sample), path("*.fastq") , emit: sim_fastq
  script:
    sample = "simulated_data" + index
    """
    art_illumina -sam -i ${reference_path}/genome.fa -p -l 150 -f 100 -m 200 -s 10 -o ${sample} 
    cp ${sample}1.fq ${sample}_R1.fastq
    cp ${sample}2.fq ${sample}_R2.fastq
    """
}

process repairFastq {
  container "$params.bbmap"
  errorStrategy 'retry'
  label 'low_mem'
  publishDir "$params.out_path/organism/fastq", mode : "copy"
  input:
    each path(ebv_dedup_reads)
  output:
    path "*_reorg.fq" , emit: repaired_fastq
  script:
    sample = ebv_dedup_reads[0].getSimpleName().replaceAll(/_R1_org/, "")
    """
    repair.sh in=${sample}_R1_org.fq in2=${sample}_R2_org.fq \
    out=${sample}_R1_reorg.fq out2=${sample}_R2_reorg.fq
    """
}

/*workflow {
  fastq_paths = Channel.fromFilePairs(params.fastq_paths, size: 4)
  out_path = Channel.from(params.out_path)
  main:
    runFastqc(fastq_paths)
    runMultiQC(runFastqc.out.raw_fastqc_out.collect{"$it"})
}*/