nextflow.enable.dsl=2

process spades {
  module "Singularity"
  container "$params.containers.sapdes"
  label 'high_mem'
  publishDir "$params.out_path/genome_assembly/", mode : "copy"
  input:
    tuple val(sample), trimmed_reads
    path reference
    val assembly_type
  output:
    path "${sample}", emit: genome_assembly
  script:
    forward = trimmed_reads[0]
    reverse = trimmed_reads[1]
    memory = "$task.memory" =~ /\d+/
    if (assembly_type == "reference_guided")
      """
      spades.py -t $task.cpus --trusted-contigs ${reference} --careful \
        -1 ${forward} -2 ${reverse} -m ${memory[0]}  -o ${sample}
      """
    else if (assembly_type == "de_novo")
      """
      spades.py -t $task.cpus --careful \
        -1 ${forward} -2 ${reverse} -m ${memory[0]}  -o ${sample}
      """
}

workflow {
  trim = Channel.fromFilePairs(params.trim, size: 2)
  out_path = Channel.from(params.out_path)

  main:
    spades(trim)
}

