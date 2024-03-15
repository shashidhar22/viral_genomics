nextflow.enable.dsl=2

process selectGenomes {
  container "$params.everse"
  errorStrategy 'retry'
  label 'mid_mem'
  cache false
  publishDir "$params.out_path/filtered_genomes", mode : "copy", pattern : "filtered/*.fasta.masked"
  input:
    path genome_assemblies
    path quast_report
    path study_metadata
  output:
    path "filtered/*", emit: filtered_genomes
    path "filtered/*.fasta.masked", emit: fasta_files
  script:
  """
  mkdir filtered
  Rscript $moduleDir/rscripts/filterGenomes.R -m 5 -s ${study_metadata} 1> stdout.txt 2> stderr.out
  """

}
process parsnp {
  container "staphb/prokka"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/parsnp", mode : "copy"
  input:
    path genome_assemblies
    path reference_genome
  output:
    path "parsnp/*", emit: parsnp_outputs
  script:
    memory = "$task.memory" =~ /\d+/
    """
    parsnp -r ${reference_genome} -d *.fasta* -o ./parsnp -p 5 -c --vcf 1> stdout.txt 2> stderr.txt
    """
}



/*process mafft{
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
*/


