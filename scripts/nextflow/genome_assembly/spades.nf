nextflow.enable.dsl=2

process spades {
  module 'SPAdes'
  label 'high_mem'
  publishDir "$params.out_path/output_data/genome_assembly/", mode : "copy"
  input:
    tuple val(sample), trim
  output:
    path "${sample}", emit: genome_assembly
  script:
    forward = trim [0]
    reverse = trim [1]
    memory = "$task.memory" =~ /\d+/
    """
spades.py -t 6 --trusted-contigs $params.ebv_1_ref --careful -1 ${forward} -2 ${reverse} -m ${memory[0]}  -o ${sample}
    """
    
}

workflow {
  trim = Channel.fromFilePairs(params.trim, size: 2)
  out_path = Channel.from(params.out_path)

  main:
    spades(trim)
}

