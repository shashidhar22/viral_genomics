nextflow.enable.dsl=2

process sam_to_bam {
  container "$params.samtools"
  label 'low_mem'
  publishDir "$params.out_path/alignment/bam/", mode : "copy"
  input:
    tuple val(sample), path(sam_path)
  output:
    tuple val(sample), path("${sample}_unmapped.bam"), emit: bam_files
  script:
    """
    samtools view -S -b -f 12 -F 256 ${sam_path} > ${sample}_unmapped.bam
    """
    
}

process bam_sort {
  container "$params.samtools"
  label 'high_mem'
  publishDir "$params.out_path/alignment/sorted/", mode : "copy", pattern: "*_sorted.ba*"
  input:
    tuple val(sample), path(bam_file)
  output:
    tuple val(sample), path("${sample}_unmapped_sorted.bam"), 
      path("${sample}_unmapped_sorted.bai"), emit: sorted_bam_file
  script:
    """
    samtools sort -n ${bam_file} -o ${sample}_unmapped_sorted.bam
    samtools index ${sample}_unmapped_sorted.bam ${sample}_unmapped_sorted.bai
    """
    
}

process idxstats {
  container "$params.samtools"
  label 'low_mem'
  publishDir "$params.out_path/idxstats/", mode : "copy"
  input:
    tuple val(sample), path(bam_file), path(index_file)
  output:
    path "${sample}_idxstats.txt", emit: idx_output
  script:
    """
    samtools idxstats ${bam_file} > ${sample}_idxstats.txt
    """
}

process coverage {
  container "$params.samtools"
  label 'low_mem'
  publishDir "$params.out_path/coverage/", mode : "copy"
  input:
    tuple val(sample), path(bam_file), path(index_file)
  output:
    path "${sample}_coverage.txt", emit: coverage_output
  script:
    """
    samtools coverage -Q 30 -q 20 -o ${sample}_coverage.txt ${bam_file}
    """
}

process keep_unaligned {
  container "$params.samtools"
  label 'low_mem'
  publishDir "$params.out_path/coverage/", mode : "copy"
  input:
    tuple val(sample), path(bam_file), path(index_file)
  output:
    tuple val(sample), path("${sample}_R*_filtered.fq.gz"), emit : filtered_reads
  script:
    """
    samtools fastq ${bam_file}\
      --threads ${task.cpus} \
      -1 ${sample}_R1_filtered.fq.gz -2 ${sample}_R2_filtered.fq.gz    
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