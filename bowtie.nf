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

// Align trimmed FASTQ reads to the host genome using Bowtie2
//
// PROCESS INPUTS:
// - Tuple with sample name, FASTQ path and path to Bowtie2 index
//  reads (Source: bbduk.out.trimmed_reads.combine(index_path))
//
// PROCESS OUTPUTS:
// - SAM file containing reads aligned to the host genome
//
// EMITTED CHANNELS:
// - sam_path :- Path to raw sam file with reads mapped to the host genome
//
// NOTE: It is better to always modularize in nextflow because of chaching 
// issues. Tuple dont get cached, incorporate sample names into the output 
// when
// TODO: 
process hostRemoval {
  container "$params.bowtie2"
  label 'mid_mem'
  errorStrategy 'retry'
  time '3d'
  publishDir "$params.out_path/alignment_reports/", mode: "copy", pattern : "*.txt"

  input:
    each fastq_path 
    path bowtie_index
    path org_index
  output:
    path "*.sam", emit: sam_path
  script:
    def r1 = fastq_path[0]
    def r2 = fastq_path[1]
    def sample = r1.getSimpleName().replaceAll(/_R1/, "")
    """
    bowtie2 -p 60 --very-fast-local -x ${bowtie_index}/genome -1 ${r1} -2 ${r2} \
      -S ${sample}.sam > stdout.log 2> ${sample}_host_alignment.txt  
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
  time '2d'
  publishDir "$params.out_path/host/alignment/", mode : "copy", pattern : "*.ba*"
  publishDir "$params.out_path/organism/unaligned_reads/", mode: "copy", pattern : "*.fq.gz"
  input:
    tuple val(sample), path(fastq_path), path(bowtie_index)
  output:
    tuple val(sample), path("${sample}_*_filtered.fq.gz"), emit: unaligned_reads
    path "${sample}_host_sorted.ba*", emit: host_bam
    path "${sample}_${ref_type}_alignment.txt", emit: alignment_stats
  script:
    r1 = fastq_path[0]
    r2 = fastq_path[1]
    """
    bowtie2 --very-sensitive-local -x ${bowtie_index}/genome -1 ${r1} -2 ${r2} \
      -S ${sample}.sam > stdout.log 2> ${sample}_${ref_type}_alignment.txt  
    samtools view -S -b -f 12 -F 256 ${sam_path} > ${sample}_unmapped.bam
    samtools view -S -b -F 1280 ${sam_path} > ${sample}_host_mapped.bam
    samtools sort  ${sample}_unmapped.bam -o ${sample}_unmapped_sorted.bam
    samtools sort  ${sample}_host_mapped.bam -o ${sample}_host_sorted.bam
    samtools index ${sample}_unmapped_sorted.bam ${sample}_unmapped_sorted.bai
    samtools index ${sample}_host_sorted.bam ${sample}_host_sorted.bai
    samtools fastq ${sample}_unmapped_sorted.bam \
      --threads ${task.cpus} \
      -1 ${sample}_R1_filtered.fq.gz -2 ${sample}_R2_filtered.fq.gz > stdout.log 2> stderr.log
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
  time '3d'
  label 'high_mem'
  publishDir "$params.out_path/alignment/ebv/bam", mode : "copy"
  input:
    each fastq_path
    path bowtie_index

  output:
    path "*sorted.ba*", emit: bam_path
  script:
    def r1 = fastq_path[0]
    def r2 = fastq_path[1]
    def sample = r1.getSimpleName()
    """
    bowtie2 -p 60 --very-sensitive-local -x ${bowtie_index}/genome -1 ${r1} -2 ${r2} \
      -S ${sample}.sam > stdout.log 2> stderr.log  
    samtools view -@ 20 -f 3 -S -b ${sample}.sam > ${sample}.bam
    samtools sort -@ 20 ${sample}.bam -o ${sample}_sorted.bam
    samtools index -@ 20 ${sample}_sorted.bam ${sample}_sorted.bai
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
  time '1d 6h'
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
