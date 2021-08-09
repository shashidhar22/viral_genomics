nextflow.enable.dsl=2

process spades {
  container "$params.spades"
  //module "SPAdes"
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.out_path/genome_assembly/$params.assembly_type", mode : "copy", pattern : "${sample}"
  publishDir "$params.out_path/contigs/$params.assembly_type/", mode : "copy", pattern : "${sample}.fasta"
  input:
    tuple val(sample), path(trimmed_reads), path(org_reference)
  output:
    tuple val(sample), path("${sample}/contigs.fasta"), path("${sample}/scaffolds.fasta"), emit: genome_assembly
    path "${sample}.fasta", emit: contigs
  script:
    forward = trimmed_reads[0]
    reverse = trimmed_reads[1]
    memory = "$task.memory" =~ /\d+/
    if ("$params.assembly_type" == "reference_guided")
      """
      spades.py -k 11,21,33,55,77 -t $task.cpus --trusted-contigs ${org_reference}/genome.fa --careful \
        -1 ${forward} -2 ${reverse} -m ${memory[0]}  -o ${sample} > stdout.log 2> stderr.log  
      cp ${sample}/contigs.fasta ${sample}.fasta
      """
    else if ("$params.assembly_type" == "de_novo")
      """
      spades.py -k 11,21,33,55,77 -t $task.cpus --careful \
        -1 ${forward} -2 ${reverse} -m ${memory[0]}  -o ${sample} > stdout.log 2> stderr.log
      cp ${sample}/contigs.fasta ${sample}.fasta
      """
}

process unicycler {
  container "$params.unicycler"
  errorStrategy 'retry'
  label 'high_mem'
  publishDir "$params.out_path/unicycler", mode : "copy", pattern : "${sample}"
  publishDir "$params.out_path/contigs/unicycler/", mode : "copy", pattern : "${sample}.fasta"
  input:
    tuple val(sample), path(trimmed_reads)
  output:
    tuple val(sample), path("${sample}/assembly.fasta"), path("${sample}/assembly.gfa"), emit: unicycler
    path "${sample}.fasta", emit: contigs
  script:
    forward = trimmed_reads[0]
    reverse = trimmed_reads[1]
    memory = "$task.memory" =~ /\d+/
    """
    unicycler -1 ${forward} -2 ${reverse} --verbosity 0 -o ${sample} > stdout.log 2> stderr.log
    cp ${sample}/assembly.fasta ${sample}.fasta
    cp ${sample}/assembly.gfa ${sample}.gfa
    """
}

process abacas {
  container "$params.abacas"
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.out_path/contigs/abacas/", mode : "copy", pattern : "${sample}.fasta"
  input:
    tuple val(sample), path(unicycler), path(unicycler_gfa), path(org_reference)
  output:
    path "${sample}.fasta", emit: contigs
  script:
    """
    abacas.pl -r ${org_reference}/genome.fa -q ${unicycler} -p nucmer > stdout.log 2> stderr.log
    awk '/^>/ {printf(">%s\\n",${sample});next;} {print}' assembly.fasta_genome.fa.fasta > ${sample}.fasta 
    """
}




workflow {
  trim = Channel.fromFilePairs(params.trim, size: 2)
  out_path = Channel.from(params.out_path)

  main:
    spades(trimmed_reads, reference, "reference_guided")
}

