nextflow.enable.dsl=2

process sam_to_bam {
  container "$params.samtools"
  errorStrategy 'retry'
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
  errorStrategy 'retry'
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
  errorStrategy 'retry'
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
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/coverage/${organism}", mode : "copy"
  input:
    tuple val(sample), path(sam_file), val(organism)
  output:
    path "${sample}_coverage.txt", emit: coverage_output
  script:
    """
    samtools view -S -b ${sam_file} > ${sample}.bam
    samtools sort ${sample}.bam -o ${sample}_sorted.bam
    samtools index ${sample}_sorted.bam ${sample}_sorted.bai
    samtools coverage -Q 30 -q 20 -o ${sample}_coverage.txt ${sample}_sorted.bam
    """
}

process keep_unaligned {
  container "$params.samtools"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/unaligned/", mode : "copy"
  input:
    tuple val(sample), path(bam_file), path(index_file)
  output:
    tuple val(sample), path("${sample}_R*_filtered.fq.gz"), emit : filtered_reads
  script:
    """
    samtools fastq ${bam_file}\
      --threads ${task.cpus} \
      -1 ${sample}_R1_filtered.fq.gz -2 ${sample}_R2_filtered.fq.gz > stdout.log 2> stderr.log
    """
}

process mpileup {
  container "$params.samtools"
  errorStrategy 'retry'
  label 'high_mem'
  publishDir "$params.out_path/mpileup/", mode : "copy"
  input:
    tuple val(sample), path(bam_file), path(index_file), path(reference)
  output:
    tuple val(sample), path("${sample}_mpileup.txt"), emit: mpileup_output
  script:
    """
    samtools mpileup -f ${reference}/genome.fa -q 20 -Q 30 -a -o ${sample}_mpileup.txt ${bam_file}
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