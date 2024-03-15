nextflow.enable.dsl=2

process prokka {
  container "$params.prokka"
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.out_path/genome_annotation/${assembly_mode}/${sample}", 
    mode : "copy"
  input:
    each path(assembled_genome)
    val assembly_mode
  output:
    path "${sample}_${assembly_mode}", emit: gene_annotation
  script:
    sample = assembled_genome.getSimpleName()
    """
    mkdir -p temp_storage
    export TMPDIR=$workDir
    prokka --outdir ${sample}_${assembly_mode} --prodigaltf $params.prodigal_trn --proteins $params.ebv_genbank --metagenome --kingdom Viruses --genus Lymphocryptovirus --rfam --prefix ${sample} --locustag EBV ${assembled_genome} > stdout.log 2> stderr.log  
    """
}

process find_epitope {
  container "$params.blast"
  label 'local_run'
  errorStrategy 'retry'
  publishDir "$params.out_path/epitope_search/", 
    mode : "copy"
  input:
    each path(protein_fasta)
    path epitope_db
  output:
    path "${sample}_blast.tsv", emit: blast_output
  script:
    sample = protein_fasta.getSimpleName()
    out_path = protein_fasta.getSimpleName()
    """
    makeblastdb -in ${protein_fasta} -dbtype nucl -out ${protein_fasta} -title ${sample}
    tblastn -db ${protein_fasta} -query ${epitope_db} -outfmt "6 sseqid qseqid evalue bitscore pident sseq qseq length qlen sstart send mismatch gaps" -out ${sample}_blast.tsv
    """
}