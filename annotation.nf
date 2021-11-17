nextflow.enable.dsl=2

process prokka {
  container "$params.prokka"
  //module "SPAdes"
  label 'high_mem'
  scratch '/fh/scratch/delete30/warren_h/ebv_enktl'
  errorStrategy 'retry'
  publishDir "$params.out_path/genome_annotation/${assembly_mode}/${sample}", 
    mode : "copy"
  input:
    each path(assembled_genome)
    val assembly_mode
  output:
    path "${sample}", emit: gene_annotation
  script:
    sample = assembled_genome.getSimpleName()
    """
    mkdir -p temp_storage
    export TMPDIR=$workDir
    prokka --outdir ${sample} --kingdom Viruses --genus Lymphocryptovirus --rfam --mincontiglen 200 --usegenus ${assembled_genome} > stdout.log 2> stderr.log  
    """
}