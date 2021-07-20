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
      path("${sample}.bai"), emit: bam_files
  script:
    sample = bowtie_align[0]
    sam_path = bowtie_align[1]
    """
    samtools view -S -b ${sam_path} > ${sample}.bam
    samtools index ${sample}.bam ${sample}.bai
    """
    
}

process bam_sort {
  module "SAMtools"
//  module "Singularity"
//  container "$params.containers.samtools"
  label 'high_mem'
  publishDir "$params.out_path/alignment/sorted/", mode : "copy", pattern: "*_sorted.ba*"
  input:
    each bam_path
  output:
    tuple val(sample), path("${sample}_sorted.bam"), 
      path("${sample}_sorted.bai"), emit: sorted_bam_file
  script:
    sample = bam_path[0]
    bam_file = bam_path[1]
    bam_outpath = bam_path[2]
    """
    samtools sort ${bam_file} -o ${sample}_sorted.bam
    samtools index ${sample}_sorted.bam ${sample}_sorted.bai
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