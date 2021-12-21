nextflow.enable.dsl=2

process spades {
  container "$params.spades"
  //module "SPAdes"
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.out_path/genome_assembly/$params.assembly_type/${sample}", 
    mode : "copy", pattern : "${sample}/*fa*"
  publishDir "$params.out_path/contigs/$params.assembly_type/", mode : "copy", 
    pattern : "${sample}.fasta"
  input:
    tuple val(sample), path(trimmed_reads)
  output:
    tuple val(sample), path("${sample}/contigs.fasta"), 
      path("${sample}/scaffolds.fasta"), path("${sample}/assembly_graph.fastg"),
      path("${sample}/assembly_graph_with_scaffolds.gfa"), emit: genome_assembly
    path "${sample}.fasta", emit: contigs
  script:
    def forward = trimmed_reads[0]
    def reverse = trimmed_reads[1]
    def memory = "$task.memory" =~ /\d+/
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
  publishDir "$params.out_path/genome_assembly/unicycler/${sample}", 
    mode : "copy", pattern : "${sample}/assembly.*"
  publishDir "$params.out_path/contigs/unicycler/", mode : "copy", 
    pattern : "${sample}.fasta"
  input:
    tuple val(sample), path(trimmed_reads)
  output:
    tuple val(sample), path("${sample}/assembly.fasta"), 
    path("${sample}/assembly.gfa"), emit: unicycler
    path "${sample}.fasta", emit: contigs
  script:
    def forward = trimmed_reads[0]
    def reverse = trimmed_reads[1]
    def memory = "$task.memory" =~ /\d+/
    """
    unicycler -1 ${forward} -2 ${reverse} --mode bold --verbosity 0 -o ${sample} > stdout.log 2> stderr.log
    cp ${sample}/assembly.fasta ${sample}.fasta
    cp ${sample}/assembly.gfa ${sample}.gfa
    """
}

process abacas {
  container "$params.abacas"
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.out_path/contigs/abacas/${mode}", mode : "copy", pattern : "${sample}_abacas.fasta"
  input:
    each path(assembly)
    val mode
    path org_reference
  output:
    path "${sample}_abacas.fasta", emit: abacas_contigs
  script:
    sample = assembly.getSimpleName()
    """
    abacas.pl -r ${org_reference}/genome.fa -q ${assembly} -p nucmer > stdout.log 2> stderr.log
    awk '/^>/ {printf(">%s\\n","${sample}");next;} {print}' ${sample}.fasta_genome.fa.fasta > ${sample}_abacas.fasta 
    """
}

process repeat_masker {
  container "$params.repeatmasker"
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.out_path/repeatmasker", mode : "copy", pattern : "${sample}.fasta*"
  publishDir "$params.out_path/contigs/repeat_masked", mode : "copy", pattern : "${sample}.fasta.masked", saveAs: { "${sample}_masked.fasta"}
  input:
    each abacas_path
    path repeat_table 
  output:
    path "${sample}.fasta.masked", emit: masked_fasta
    path "${sample}_masked.fasta", emit: contigs
    path "${sample}.fasta.cat", emit: repeat_catalog
    path "${sample}.fasta.out", emit: repeat_output
    path "${sample}.fasta.tbl", emit: repeat_table
  script:
    sample = abacas_path.getBaseName()
    """
    RepeatMasker -pa 4 -lib ${repeat_table} -s -e rmblast -dir ./ ${abacas_path} > stdout.log 2> stderr.log
    cp ${sample}.fasta.masked ${sample}_masked.fasta
    """
}

workflow {
  trim = Channel.fromFilePairs(params.trim, size: 2)
  out_path = Channel.from(params.out_path)

  main:
    spades(trimmed_reads, reference, "reference_guided")
}

