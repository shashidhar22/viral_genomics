nextflow.enable.dsl=2
process bbduk {
  //module "Singularity"
  //container "$params.containers.bbmap"
  module "BBMap"
  label 'mid_mem'
  publishDir "$params.out_path/trimmed_reads/", mode : "copy"
  input:
    each fastq_path
    path adapters
  output:
    tuple val(sample), path("*clean.fq"), emit: trimmed_reads
  script:
    sample = fastq_path[0]
    r1l1 = fastq_path[1][0]
    r2l1 = fastq_path[1][1]
    r1l2 = fastq_path[1][2]
    r2l2 = fastq_path[1][3]
    """
    zcat ${r1l1} ${r1l2} > ${sample}_R1.fastq 
    zcat ${r2l1} ${r2l2} > ${sample}_R2.fastq 
    bbduk.sh in1=${sample}_R1.fastq in2=${sample}_R2.fastq \
    out1=${sample}_R1.clean.fq out2=${sample}_R2.clean.fq \
    ref=${adapters}\
    ktrim=r k=23 mink=11 hdist=1 tpe=f tbo=t
    """
    
}

/*workflow {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: -1)
  out_path = Channel.from(params.out_path)
  no_human = Channel.fromPath(params.no_human)

  main:
    bbduk(fastq_path)
}
*/
