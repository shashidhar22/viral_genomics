nextflow.enable.dsl=2

process biosino_download {
  label 'high_mem'
  errorStrategy 'retry'
  input:
    each row
  output:
    tuple val(sample), path("${sample}_*.fq.gz"), emit: fastq_files
  script:
    def sample = row.sample_id + '_' + row.fileName.replaceAll("_{1,2}.fq.gz", "")
    def read_name = row.fileName
    def read_url = row.url
    def read_out = row.sample_id + '_' + row.fileName
    """
    wget -c --retry-connrefused --retry-on-http-error=552,502 -t 0 -T 600 -O ${read_out} ${read_url}
    """
}