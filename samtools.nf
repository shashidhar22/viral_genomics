nextflow.enable.dsl=2


// Extract host aligned and organism aligned reads from the raw SAM file 
// and store it as BAM files
//
// PROCESS INPUTS:
// - Raw SAM file (hostRemoval.out.sam_path)
//
// PROCESS OUTPUTS:
// - BAM file containing Host aligned reads with index
// - BAM file containing reads that didnt map to the Host with index
//
// EMITTED CHANNELS:
// - host_bam_path :- List containing host BAM and BAI files
// - unmp_bam_path :- List containing unmapped reads BAM and BAI files
//
// NOTE: 
//
// TODO: 

process splitSam {
  container "$params.samtools"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/host/alignments", mode : "copy",
    pattern : "*_host.ba*"
  publishDir "$params.out_path/unmapped/alignments", mode : "copy",
    pattern : "*_unmp.ba*"
  input:
    each path(sam_path)
  output:
    path "*_host.ba*", emit : host_bam_path
    path "*_unmp.ba*", emit : org_bam_path
  script:
    def sample = sam_path.getSimpleName()
    """
    samtools view -@ 20 -S -b -f 12 -F 256 ${sample}.sam >> ${sample}_unmapped.bam
    samtools view -@ 20 -S -b -F 1280 ${sample}.sam >> ${sample}_host_reads.bam
    samtools sort -@ 20 ${sample}_unmapped.bam -o ${sample}_unmp.bam
    samtools sort -@ 20 ${sample}_host_reads.bam -o ${sample}_host.bam
    samtools index -@ 20 ${sample}_unmp.bam ${sample}_unmp.bai
    samtools index -@ 20 ${sample}_host.bam ${sample}_host.bai
    """
    
}

process extractHostReads {
  container "$params.samtools"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/host/fastq", mode : "copy",
    pattern : "*_host.fq"
  input:
    path sam_path
  output:
    path "*_host.fq", emit : host_fastq_files
  script:
    def sample = sam_path[0].getSimpleName().replaceAll(/_host/, "")
    """
    samtools fastq -f 3 ${sample}_host.bam \
      --threads ${task.cpus} \
      -1 ${sample}_R1_host.fq -2 ${sample}_R2_host.fq > stdout.log 2> stderr.log    
    """
}

process extractOrgReads {
  container "$params.bowtie2"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/organism/fastq", mode : "copy",
    pattern : "*_org.fq"
  publishDir "$params.out_path/organism/alignments", mode : "copy",
    pattern : "*_org.ba*"
  publishDir "$params.out_path/unmapped/fastq", mode : "copy",
    pattern : "*_unmp.fq"
  input:
    each path(sam_path)
    path org_index
  output:
    path "*_org.fq", emit : org_fastq_files
    tuple val(sample), path("*_org.fq"), emit : org_fastqc_input
    path "*_org.ba*" , emit: org_bam_path
  script:
    sample = sam_path[0].getSimpleName().replaceAll(/_unmp/, "")
    """
    samtools fastq -G 1024 ${sample}_unmp.bam \
      --threads ${task.cpus} \
      -1 ${sample}_R1_unmp.fq -2 ${sample}_R2_unmp.fq > stdout.log 2> stderr.log    
    bowtie2 -p 20 --very-sensitive-local -x ${org_index}/genome \
      -1 ${sample}_R1_unmp.fq -2 ${sample}_R2_unmp.fq \
      -S ${sample}_org_mapped.sam >> stdout.log 2>> ${sample}_org_alignment.txt
    samtools view -@ 20 -S -b ${sample}_org_mapped.sam > ${sample}_org_mapped.bam
    samtools fixmate -@ 20 -m ${sample}_org_mapped.bam ${sample}_org_fixmate.bam
    samtools sort -@ 20 -o ${sample}_org_sorted.bam ${sample}_org_fixmate.bam
    samtools markdup -@ 20 -r ${sample}_org_sorted.bam ${sample}_org.bam
    samtools index -@ 20 ${sample}_org.bam ${sample}_org.bai
    samtools fastq -f 3 ${sample}_org.bam \
      --threads ${task.cpus} \
      -1 ${sample}_R1_org.fq -2 ${sample}_R2_org.fq \
      >> stdout.log 2>> stderr.log
    
    """
}



process bam_sort {
  container "$params.samtools"
  errorStrategy 'retry'
  label 'high_mem'
  time '1d 6h'
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
  time '1d 6h'
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
  time '1d 6h'
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
  time '1d 6h'
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
  time '1d 6h'
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