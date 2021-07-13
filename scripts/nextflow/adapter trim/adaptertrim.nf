nextflow.enable.dsl=2

process bbduk {
  module 'BBMap'
  label 'low_mem'
  publishDir "$params.out_path/input_data/trim_human/", mode : "copy"
  input:
    tuple val(sample), path(fastq_path)
  output:
    path "*clean.fq", emit: trim
  script:
    r1l1 = fastq_path[0]
    r2l1 = fastq_path[1]
    r1l2 = fastq_path[2]
    r2l2 = fastq_path[3]
    """
    zcat ${r1l1} ${r1l2} > ${sample}_R1.fastq 
    zcat ${r2l1} ${r2l2} > ${sample}_R2.fastq 
    bbduk.sh in1=${sample}_R1.fastq in2=${sample}_R2.fastq out1=${sample}_R1.clean.fq out2=${sample}_R2.clean.fq ref=/fh/fast/warren_h/users/ebrierch/ebv_seq_project/input_data/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
    """
    
}

workflow {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: -1)
  out_path = Channel.from(params.out_path)
  no_human = Channel.fromPath(params.no_human)

  main:
    bbduk(fastq_path)
}
