nextflow.enable.dsl=2
include { bbduk } from './trim_reads'
include { runFastqc; runMultiQC } from './dataqc'
include { bowtie_Align } from './bowtie'
include { sam_to_bam;bam_sort;coverage } from './samview'

workflow alignAndCoverage {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: 4)
  adapter_path = Channel.fromPath(params.adapters)
  out_path = Channel.fromPath(params.out_path)
  index_path = Channel.fromPath(params.index.host.bowtie.genome)
  reference_path = Channel.from(params.references.host.genome)
  main:
    bbduk(fastq_path, adapter_path)
    bowtie_Align(bbduk.out.trimmed_reads, index_path)
    sam_to_bam(bowtie_Align.out.sam_path)
    bam_sort(sam_to_bam.out.bam_files)
    coverage(bam_sort.out.sorted_bam_file)
}

workflow trimReadsandQC {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: 4)
  adapter_path = Channel.fromPath(params.adapters)

  out_path = Channel.from(params.out_path)

  main:
    bbduk(fastq_path, adapter_path)
    runFastqc(bbduk.out.trimmed_reads, "trimmed")
    runMultiQC(runFastqc.out.fastqc_results.collect{"$it"})

}

workflow rawFastqQC {
  fastq_paths = Channel.fromFilePairs(params.fastq_paths, size: 4)
  out_path = Channel.from(params.out_path)

  main:
    runFastqc(fastq_paths, "raw")
    runMultiQC(runFastqc.out.fastqc_results.collect{"$it"})
}

