nextflow.enable.dsl=2

process spades {
  module 'SPAdes'
  label 'high_mem'
  publishDir "$params.out_path/output_data/genome_assembly_no_ref/", mode : "copy"
  input:
    tuple val(sample), path(trim)
  output:
    path "${sample}", emit: genome_assembly
  script:
    forward = trim [0]
    reverse = trim [1]
    memory = "$task.memory" =~ /\d+/
    """
spades.py -t 6 --careful -1 ${forward} -2 ${reverse} -m ${memory[0]}  -o ${sample}
    """
    
}

process assembly_QC {
  module 'QUAST'
  label 'low_mem'
  publishDir "$params.out_path/output_data/genome_assembly_no_ref/quast", mode : "copy"
  input:
    each genome_assembly
  output:
    path "${sample}", emit: quast
  script:
    memory = "$task.memory" =~ /\d+/
    sample = genome_assembly.getParent().getSimpleName()
    """
    quast.py ${genome_assembly} -r $params.ebv_1_ref  -o ${sample}
    """
    
}

process runMultiQC {
  module 'MultiQC'
  label 'low_mem'
  publishDir "$params.out_path/output_data/genome_assembly_no_ref/multiqc", mode : "copy"
  input:
    path quast
  output:
    path "multiqc_report.html", emit: mutli_qc_out
  script:
    """
    multiqc .
    """
    
}
workflow {
  trim = Channel.fromFilePairs(params.trim, size: -1)
  out_path = Channel.from(params.out_path)

  main:
    spades(trim)
    assembly_QC(spades.out.genome_assembly)
    runMultiQC(assembly_QC.out.quast)
}