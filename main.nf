nextflow.enable.dsl=2
include { bbduk ; bbnorm } from './trim_reads'
include { runFastqc; runMultiQC; quast; quast as quast_masker; quast as quast_uniclyer; quast as quast_abacas ; getFastq ; generateFastq} from './dataqc'
include { bowtie_Align; align_ebv; align_filtered_ebv ; align_filtered_ebv as align_sampled_ebv} from './bowtie'
include { sam_to_bam;bam_sort;coverage;keep_unaligned;mpileup } from './samview'
include { kmer_count } from './khmer'
include { spades; unicycler; abacas; repeat_masker } from './spades'
include { coverage as ebv_coverage } from './samview'
include { haplotype_caller; cnn_score_variants ; annotate_vcfs} from './variant_calling'
include { call_host_variants; generate_gvcf_table; consolidate_gvcfs; genotype_gvcfs ; variant_scoring } from './variant_calling'
include { prokka } from './annotation'

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

//
workflow removeHostAndAssemble {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: 4, followLinks: false)
  //public_fastq = Channel.fromSRA(params.public_list, apiKey: params.api_key)
  public_list = Channel.from(params.public_list)
  public_contigs = Channel.fromPath(params.ebv_public_data)
  ebv_repeat = Channel.fromPath(params.ebv_repeats)
  //public_fastq = Channel.fromFilePairs(params.public_fastq, size: 2, followLinks: false)
  adapter_path = Channel.fromPath(params.adapters)
  simulated_reads = Channel.from(1,2,3,4,5,6,7,8,9,10)
  out_path = Channel.fromPath(params.out_path)
  index_path = Channel.fromPath(params.index.host.bowtie.sans_unplacedebv)
  host_reference_path = Channel.from(params.references.host.sans_unplacedebv)
  org_reference_path = Channel.from(params.references.organism.ebv_1_ref)
  org_index_path = Channel.fromPath(params.index.organism.bowtie.ebv_1_ref)
  
  main:
    getFastq(public_list)
    generateFastq(simulated_reads, org_reference_path)
    bbduk(fastq_path.concat(generateFastq.out.sim_fastq).combine(adapter_path))
    bowtie_Align(bbduk.out.trimmed_reads.combine(index_path))
    coverage(bowtie_Align.out.sam_path)
    sam_to_bam(bowtie_Align.out.sam_path)
    bam_sort(sam_to_bam.out.bam_files)
    keep_unaligned(bam_sort.out.sorted_bam_file)
    align_ebv(bbduk.out.trimmed_reads.concat(generateFastq.out.sim_fastq).combine(org_index_path))
    align_filtered_ebv(keep_unaligned.out.filtered_reads.combine(org_index_path))
    haplotype_caller(align_filtered_ebv.out.sam_path.combine(org_reference_path))
    cnn_score_variants(haplotype_caller.out.vcf_paths.combine(org_reference_path))
    annotate_vcfs(cnn_score_variants.out.annotated_vcf_paths)
    bbnorm(keep_unaligned.out.filtered_reads)
    align_sampled_ebv(bbnorm.out.sampled_reads.combine(org_index_path))
    //ebv_coverage(align_filtered_ebv.out.sam_path)
    spades(bbnorm.out.sampled_reads.combine(org_reference_path))
    //quast_prep(spades.out.genome_assembly)
    quast("Spades", spades.out.contigs.toSortedList(), org_reference_path)
    unicycler(seqtk.out.sampled_reads)
    quast_uniclyer("Unicycler", unicycler.out.contigs.toSortedList().unique(), org_reference_path)
    abacas(unicycler.out.unicycler.combine(org_reference_path))
    quast_abacas("Abacas", abacas.out.contigs.toSortedList(), org_reference_path)    
    repeat_masker(abacas.out.contigs.concat(public_contigs), ebv_repeat)
    quast_masker("RepeatMasker", repeat_masker.out.contigs.toSortedList(), org_reference_path)
    mpileup(align_sampled_ebv.out.sam_path.combine(org_reference_path))
    prokka(unicycler.out.contigs)
}


workflow hostVariantCalling {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: 4, followLinks: false)
  adapter_path = Channel.fromPath(params.adapters)
  out_path = Channel.fromPath(params.out_path)
  index_path = Channel.fromPath(params.index.host.bowtie.sans_unplacedebv)
  host_reference_path = Channel.from(params.references.host.sans_unplacedebv)
  varcall_host_path = Channel.from(params.references.host.variant_calling)

  main:
    bbduk(fastq_path.combine(adapter_path))
    bowtie_Align(bbduk.out.trimmed_reads.combine(index_path))
    coverage(bowtie_Align.out.sam_path)
    sam_to_bam(bowtie_Align.out.sam_path)
    bam_sort(sam_to_bam.out.bam_files)
    call_host_variants(bam_sort.out.sorted_bam_file.combine(varcall_host_path))
    generate_gvcf_table(call_host_variants.out.gvcf_paths, call_host_variants.out.gvcf_index_paths)
    consolidate_gvcfs(call_host_variants.out.gvcf_paths, call_host_variants.out.gvcf_index_paths, generate_gvcf_table.out.gvcf_map_file)
    genotype_gvcfs(consolidate_gvcfs.out.gvcf_database, varcall_host_path)
    variant_scoring(genotype_gvcfs.out.joint_vcf, genotype_gvcfs.out.joint_vcf_index, varcall_host_path)
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