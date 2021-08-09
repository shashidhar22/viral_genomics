nextflow.enable.dsl=2
include { bbduk } from './trim_reads'
include { runFastqc; runMultiQC; quast; quast as quast_uniclyer; quast as quast_abacas ; getFastq} from './dataqc'
include { bowtie_Align; align_ebv} from './bowtie'
include { sam_to_bam;bam_sort;coverage;keep_unaligned } from './samview'
include { kmer_count } from './khmer'
include { spades; unicycler; abacas } from './spades'
include { coverage as ebv_coverage } from './samview'
workflow alignAndCoverage {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: 4)
  adapter_path = Channel.fromPath(params.adapters)
  out_path = Channel.fromPath(params.out_path)
  index_path = Channel.fromPath(params.index.host.bowtie.genome)
  reference_path = Channel.from(params.references.host.genome)
  main:
    bbduk(fastq_path, adapter_path)
    bowtie_Align(bbduk.out.trimmed_reads)
    sam_to_bam(bowtie_Align.out.sam_path)
    bam_sort(sam_to_bam.out.bam_files)
    coverage(bam_sort.out.sorted_bam_file)
    kmer_count(bbduk.out.trimmed_reads)
}


workflow removeHostAndAssemble {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: 4, followLinks: false)
  //public_fastq = Channel.fromSRA(params.public_list, apiKey: params.api_key)
  public_list = Channel.from(params.public_list)
  //public_fastq = Channel.fromFilePairs(params.public_fastq, size: 2, followLinks: false)
  adapter_path = Channel.fromPath(params.adapters)
  out_path = Channel.fromPath(params.out_path)
  index_path = Channel.fromPath(params.index.host.bowtie.sans_unplacedebv)
  host_reference_path = Channel.from(params.references.host.sans_unplacedebv)
  org_reference_path = Channel.from(params.references.organism.ebv_1_ref)
  org_index_path = Channel.fromPath(params.index.organism.bowtie.ebv_1_ref)
  main:
    getFastq(public_list)
    bbduk(fastq_path.combine(adapter_path))
    bowtie_Align(bbduk.out.trimmed_reads.concat(getFastq.out.public_fastq).combine(index_path))
    coverage(bowtie_Align.out.sam_path)
    sam_to_bam(bowtie_Align.out.sam_path)
    bam_sort(sam_to_bam.out.bam_files)
    keep_unaligned(bam_sort.out.sorted_bam_file)
    align_ebv(keep_unaligned.out.filtered_reads.combine(org_index_path))
    ebv_coverage(align_ebv.out.sam_path)
    spades(keep_unaligned.out.filtered_reads.combine(org_reference_path))
    //quast_prep(spades.out.genome_assembly)
    quast("Spades", spades.out.contigs.collect{"$it"}, org_reference_path)
    unicycler(keep_unaligned.out.filtered_reads)
    quast_uniclyer("Unicycler", unicycler.out.contigs.collect{"$it"}, org_reference_path)
    abacas(unicycler.out.unicycler.combine(org_reference_path))
    quast_abacas("Abacas", abacas.out.contigs.collect{"$it"}, org_reference_path)

}

workflow testSRA {
  public_fastq = Channel.fromSRA('SRR9641815', apiKey: params.api_key)
  main:
    public_fastq.view()
}

workflow trimReadsandQC {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: 4)
  adapter_path = Channel.fromPath(params.adapters, )

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

workflow spadesRun {
  fastq_paths = Channel.fromFilePairs(params.fastq_paths, size: 2)
}