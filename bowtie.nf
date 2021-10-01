nextflow.enable.dsl=2
// Create a bowtie index give a reference genome FASTA file
//
// PROCESS INPUTS:
// - Path to refrence FASTA file
//  reads (Source: reference_path)
//
// PROCESS OUTPUTS:
// - Bowtie2 index files 
// EMITTED CHANNELS:
// - bowtie_index :- Channel containing Bowtie2 index files for the reference
//  genome
//
// NOTE: 
//
// TODO: 
process indexBowtie{
  container "$params.bowtie2"
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.out_path/index/bowtie", mode : "copy"
  input:
    path reference_fasta
  output:
    path "${index_name}.*", emit: bowtie_index
  script:
    index_name = reference_fasta.getSimpleName()
    """
    bowtie2-build ${reference_fasta} ${index_name}
    """
}

// Align trimmed FASTQ reads to the reference genome using Bowtie2
//
// PROCESS INPUTS:
// - Tuple with sample name, FASTQ path and path to Bowtie2 index
//  reads (Source: bbduk.out.trimmed_reads.combine(index_path))
//
// PROCESS OUTPUTS:
// - SAM file containing reads aligned to the reference genome
// EMITTED CHANNELS:
// - sam_path :- Tuple with sample name and the path to the SAM file
//
// NOTE: 
//
// TODO: 
process bowtie_Align {
  container "$params.bowtie2"
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.out_path/alignment/sam", mode : "copy"
  input:
    tuple val(sample), path(fastq_path), path(bowtie_index)

  output:
    tuple val(sample), path("${sample}.sam"), emit: sam_path
  script:
    r1 = fastq_path[0]
    r2 = fastq_path[1]
    """
    bowtie2 --very-sensitive-local -x ${bowtie_index}/genome -1 ${r1} -2 ${r2} \
      -S ${sample}.sam > stdout.log 2> stderr.log  
    """
}

// Align trimmed FASTQ reads to the EBV genome
//
// PROCESS INPUTS:
// - Tuple with sample name, FASTQ path and path to Bowtie2 index for the 
// EBV genome
//  (Source: bbduk.out.trimmed_reads.combine(org_index_path))
//
// PROCESS OUTPUTS:
// - BAM file and index with reads mapped to the EBV genome
// EMITTED CHANNELS:
// - sam_path :- Tuple with sample name, path to BAM file and path to the 
// index file
//
// NOTE: 
//
// TODO: 
process align_ebv {
  container "$params.samtools_bowtie"
  errorStrategy 'retry'
  label 'high_mem'
  publishDir "$params.out_path/alignment/ebv/bam", mode : "copy"
  input:
    tuple val(sample), path(fastq_path), path(bowtie_index)

  output:
    tuple val(sample), path("${sample}_sorted.bam"), 
      path("${sample}_sorted.bai"), emit: sam_path
  script:
    r1 = fastq_path[0]
    r2 = fastq_path[1]
    """
    bowtie2 --very-sensitive-local -x ${bowtie_index}/genome -1 ${r1} -2 ${r2} \
      -S ${sample}.sam > stdout.log 2> stderr.log  
    samtools view -S -b ${sample}.sam > ${sample}.bam
    samtools sort ${sample}.bam -o ${sample}_sorted.bam
    samtools index ${sample}_sorted.bam ${sample}_sorted.bai
    """
}

// Align trimmed and host filtered FASTQ reads to the EBV genome
//
// PROCESS INPUTS:
// - Tuple with sample name, FASTQ path and path to Bowtie2 index for the 
// EBV genome
//  (Source: bbduk.out.trimmed_reads.combine(org_index_path))
//
// PROCESS OUTPUTS:
// - BAM file and index with reads mapped to the EBV genome
// EMITTED CHANNELS:
// - sam_path :- Tuple with sample name, path to BAM file and path to the 
// index file
//
// NOTE: 
//
// TODO: 
process align_filtered_ebv {
  container "$params.samtools_bowtie"
  errorStrategy 'retry'
  label 'high_mem'
  publishDir "$params.out_path/alignment/ebv/bam_filtered", mode : "copy"
  input:
    tuple val(sample), path(fastq_path), path(bowtie_index)

  output:
    tuple val(sample), path("${sample}_filtered_sorted.bam"), 
      path("${sample}_filtered_sorted.bai"), emit: sam_path
  script:
    r1 = fastq_path[0]
    r2 = fastq_path[1]
    """
    bowtie2 --very-sensitive-local -x ${bowtie_index}/genome -1 ${r1} -2 ${r2} \
      -S ${sample}_filtered.sam > stdout.log 2> stderr.log  
    samtools view -S -b ${sample}_filtered.sam > ${sample}_filtered.bam
    samtools sort ${sample}_filtered.bam -o ${sample}_filtered_sorted.bam
    samtools index ${sample}_filtered_sorted.bam ${sample}_filtered_sorted.bai
    """
}

workflow {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: -1)
  out_path = Channel.from(params.out_path)
  bowtie_index = Channel.from(params.bowtie_index)
  
  main:
    bowtie_Align(fastq_path)

}

