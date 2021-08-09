nextflow.enable.dsl=2

process kmer_count {
  container "$params.khmer"
  label 'high_mem'
  publishDir "$params.out_path/kmer_counts/abundace/", mode : "copy"
  input:
    tuple val(sample), path(fastq_path)
  output:
    path "${sample}_r*_abundace.tsv", emit: kmer_abundace
  script:
    r1 = fastq_path[0]
    r2 = fastq_path[1]
    """
    load-into-counting.py -x 2e8 -k 20 ${sample}.kh ${r1} 
    abundance-dist.py -s ${sample}.kh ${r1} ${sample}_r1_abundace.tsv
    load-into-counting.py -x 2e8 -k 20 ${sample}.kh ${r2} 
    abundance-dist.py -s ${sample}.kh ${r2} ${sample}_r2_abundace.tsv
    """
}