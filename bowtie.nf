nextflow.enable.dsl=2

process indexBowtie{
  container "$params.bowtie2"
  label 'high_mem'
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


process bowtie_Align {
  module "Bowtie2"
  label 'high_mem'
  publishDir "$params.out_path/alignment/sam", mode : "copy"
  input:
    tuple val(sample), path(fastq_path), path(bowtie_index)

  output:
    tuple val(sample), path("${sample}.sam"), val("NA"), emit: sam_path
  script:
    r1 = fastq_path[0]
    r2 = fastq_path[1]
    """
    bowtie2 -x ${bowtie_index}/genome -1 ${r1} -2 ${r2} -S ${sample}.sam > stdout.log 2> stderr.log  
    """
}

workflow {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: -1)
  out_path = Channel.from(params.out_path)
  bowtie_index = Channel.from(params.bowtie_index)
  
  main:
    bowtie_Align(fastq_path)

}

