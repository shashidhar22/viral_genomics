nextflow.enable.dsl=2

process spades {
  //container "$params.spades"
  module "SPAdes"
  label 'high_mem'
  publishDir "$params.out_path/genome_assembly/", mode : "copy"
  input:
    tuple val(sample), path(trimmed_reads), path(org_reference)
  output:
    path "${sample}", emit: genome_assembly
  script:
    forward = trimmed_reads[0]
    reverse = trimmed_reads[1]
    memory = "$task.memory" =~ /\d+/
    if ("$params.assembly_type" == "reference_guided")
      """
      spades.py -t $task.cpus --trusted-contigs ${org_reference} --careful \
        -1 ${forward} -2 ${reverse} -m ${memory[0]}  -o ${sample} > stdout.log 2> stderr.log  
      """
    else if ("$params.assembly_type" == "de_novo")
      """
      spades.py -t $task.cpus --careful \
        -1 ${forward} -2 ${reverse} -m ${memory[0]}  -o ${sample} > stdout.log 2> stderr.log
      """
}

workflow {
  trim = Channel.fromFilePairs(params.trim, size: 2)
  out_path = Channel.from(params.out_path)

  main:
    spades(trimmed_reads, reference, "reference_guided")
}

