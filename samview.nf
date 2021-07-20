nextflow.enable.dsl=2

process sam_to_bam {
  module "SAMtools"
//  module "Singularity"
//  container "$params.containers.samtools"
  label 'low_mem'
  publishDir "$params.out_path/alignment/bam/", mode : "copy"
  input:
    each bowtie_align
  output:
    tuple val(sample), path("${sample}.bam"), 
      val("$params.out_path/alignment/bam/"), emit: bam_files
  script:
    sample = bowtie_align[0]
    sam_path = bowtie_align[1]
    """
    samtools view -S -b ${sam_path} > ${sample}.bam
    """
    
}

process bam_sort {
  module "SAMtools"
//  module "Singularity"
//  container "$params.containers.samtools"
  label 'high_mem'
  publishDir "$params.out_path/alignment/sorted/", mode : "copy", pattern: "*.sorted.bam"
  input:
    each bam_path
  output:
    tuple val(sample), path("${sample}.sorted.bam"), 
      val("$params.out_path/alignment/sorted/"), emit: sorted_bam_file
  script:
    sample = bam_path[0]
    bam_file = bam_path[1]
    bam_outpath = bam_path[2]
    """
    samtools sort ${bam_file} -o ${sample}_sorted.bam
    """
    
}

process bam_index {
  module "SAMtools"
//  module "Singularity"
//  container "$params.containers.samtools"
  label 'high_mem'
  publishDir "${out_path}", mode : "copy", pattern: "*.bai"
  input:
    each bam_path
  output:
    tuple val(sample), path(bam_file), path("${sample}.bai"), emit: indexed_path

  script:
    sample = bam_path[0]
    bam_file = bam_path[1]
    out_path = bam_path[2]
    """
    samtools index ${bam_file} ${sample}.bai
    """
}

process keep_unaligned {
  module "SAMtools"
//  module "Singularity"
//  container "$params.containers.samtools"
  label 'low_mem'
  publishDir "$params.out_path/unaligned_reads/", mode : "copy"
  input:
    each bam_path
  output:
    path "${sample}_R*.filtered.fq.gz", emit: unaligned_fastq
  script:
    sample = bam_path[0]
    bam_file = bam_path[1]
    index_file = bam_path[2]
    """
    samtools fastq ${bam_file} --threads 8 -F 12 -1 ${sample}_R1_filtered.fq.gz \
    -2 ${sample}_R2_filtered.fq.gz
    """
    
}


process idxstats {
  module "SAMtools"
//  module "Singularity"
//  container "$params.containers.samtools"
  label 'low_mem'
  publishDir "$params.out_path/idxstats/", mode : "copy"
  input:
    each bam_path
  output:
    path "${sample}_idxstats.txt", emit: idx_output
  script:
    sample = bam_path[0]
    bam_file = bam_path[1]
    out_path = bam_path[2]
    """
    samtools idxstats ${bam_file} > ${sample}_idxstats.txt
    """
}

process coverage {
  module "SAMtools"
//  module "Singularity"
//  container "$params.containers.samtools"
  label 'low_mem'
  publishDir "$params.out_path/coverage/", mode : "copy"
  input:
    each bam_path
  output:
    path "${sample}_coverage.txt", emit: coverage_output
  script:
    sample = bam_path[0]
    bam_file = bam_path[1]
    index_path = bam_path[2]
    """
    samtools coverage -Q 30 -q 20 -o ${sample}_coverage.txt ${bam_file}
    """
}

workflow {
  bowtie_align = Channel.fromPath(params.bowtie_align)
  out_path = Channel.from(params.out_path)

  main:
    sam_bam(bowtie_align)
    keep_unaligned(bowtie_align)
    bam_sort(sam_bam.out.bam_files)
    sam_index(bam_sort.out.sorted_bam_file)
    idxstats(bam_sort.out.sorted_bam_file)
}