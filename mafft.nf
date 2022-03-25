nextflow.enable.dsl=2




process mafft{
  module "Singularity"
  container "$params.container.maaft"
  label 'high_mem'
  publishDir "$params.out_path/mafft", mode : "copy"
  input:
    path assembled_genomes
  output:
    path "mafft_alignment.mafft", emit: mafft
  script:
    
    """
    zcat 
    mafft mafft_alignment.fasta > mafft_alignment.mafft
    """
}

workflow gblocks {
  module "Singularity"
  container "params.container.glbocks"
  label 'low_mem'

  publishDir "$params.out_path/gblock", mode: 'copy'

  input:
    file mafft_alignment

  output:
    file "*-gb", emit: glbocks

  script:
    """
    Gblocks ${mafft_alignment} -t=d -p=n -b3=8 -b4=10 -b5=h
    """
}

workflow {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: -1)
  out_path = Channel.from(params.out_path)
  bowtie_index = Channel.from(params.bowtie_index)
  
  main:
    bowtie_Align(fastq_path)

}



mafft /fh/fast/warren_h/ebv_enktl/results/sans_unplaced_ebv_assembly_and_annotation/contigs/cat_repeat_masked_limited.fasta > /fh/fast/warren_h/ebv_enktl/results/sans_unplaced_ebv_assembly_and_annotation/contigs/mafft/repeat_masked_limited.fasta.mafft