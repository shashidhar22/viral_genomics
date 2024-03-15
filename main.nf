nextflow.enable.dsl=2
include { runFastQC as fastqQCPublic ; runFastQC as fastqQCRaw ; runFastQC as fastqQCTrim ; runFastQC as fastQCEbv} from './dataqc'
include { bbduk ; bbnorm ; bbduk as bbdukWGS } from './trim_reads'
include { hostRemoval ; align_ebv as alignEBV } from './bowtie'
include { repairFastq } from './dataqc'
include { abacas as abacasSpades ; abacas as abacasUnicycler} from './spades'
include { prokka as prokkaSpades ; prokka as prokkaUnicycler ; prokka as prokkaCombined } from './annotation'
include { quast as quastSpades ; quast as quastUnicycler ; quast as quastCombined } from './dataqc'
include { haplotype_caller as haplotypeCaller ; haplotype_caller_wgs as haplotypeCallerWGS} from './variant_calling'
include { cnn_score_variants as cnnScoreVariants ; cnn_score_variants as cnnScoreVariantsWGS} from './variant_calling'
include { annotate_vcfs as annotateVCFs ; annotate_vcfs as annotateVCFsWGS} from './variant_calling'
include { biosino_download } from './data_download'
include { runFastQC; runMultiQC; quast; quast as quast_masker; quast as quast_uniclyer; quast as quast_abacas ; getFastq ; generateFastq} from './dataqc'
//include { bowtie_Align; align_ebv; align_filtered_ebv ; align_filtered_ebv as align_sampled_ebv} from './bowtie'
include { bam_sort;coverage;keep_unaligned;mpileup } from './samtools'
include { spades; unicycler; abacas; repeat_masker } from './spades'
include { coverage as ebv_coverage } from './samtools'
include { splitSam ; extractHostReads ; extractOrgReads } from './samtools'
include { create_alternate_reference as createAlternateReference ; create_alternate_reference_wgs as createAlternateReferenceWGS} from './variant_calling'
include { call_host_variants; generate_gvcf_table; consolidate_gvcfs; genotype_gvcfs ; variant_scoring ; variant_call_muscle } from './variant_calling'
include { prokka ; prokka as prokka_spades ; prokka as prokka_abacas ; find_epitope} from './annotation'
include { parsnp ; selectGenomes } from './phylogeny'
workflow ebvAssembly {
  ed_data = Channel.fromFilePairs(params.ed_data, size: 4)
  elshafa_data = Channel.fromFilePairs(params.elshafa_data, size: 2)
  fastq_path = ed_data.concat(elshafa_data)
  public_list = Channel.from(params.public_list)
  adapter_path = Channel.fromPath(params.adapters)
  public_fasta = Channel.fromPath(params.assembled_genome )
  out_path = Channel.fromPath(params.out_path)
  ebv_reference = Channel.fromPath(params.references.organism.ebv_1_ref + '/genome.fa')
  ebv_masked_reference = Channel.fromPath(params.references.masked.ebv_1_ref + '/genome.fa')
  human_reference = Channel.fromPath(params.references.host.sans_ebv + '/genome.fa')
  ebv_index = Channel.fromPath(params.references.organism.ebv_1_ref)
  ebv_masked_index = Channel.fromPath(params.references.masked.ebv_1_ref)
  human_index = Channel.fromPath(params.references.host.sans_ebv)
  multi_config = Channel.fromPath(params.multiqc_config)
  duplex_vcf = Channel.fromPath(params.duplex_vcfs)
  duplex_tabix = Channel.fromPath(params.duplex_tabix)
  duplex_list = Channel.from(params.duplex_list)
  ebv_findex = Channel.fromPath(params.references.masked.ebv_1_ref + '/genome.fa.fai')
  ebv_dict = Channel.fromPath(params.references.masked.ebv_1_ref + '/genome.dict')
  out_path = Channel.fromPath(params.out_path)
  repeat_table = Channel.fromPath(params.ebv_repeats)
  wgs_path = Channel.fromFilePairs(params.biosino_data,  size: 2)
  epitope_db = Channel.fromPath(params.ebv_epitope_db)
  main:
    getFastq(public_list)
    //Run FastQC on raw read
    fastqQCRaw(fastq_path.concat(getFastq.out.public_fastq), "public")
    // Download biosino data
    //biosino_download(biosino_data.splitCsv(header: true, sep: '\t'))
    //fastqQCPublic(biosino_download.out.fastq_files.groupTuple(), "public")
    // Merge fastq channels
    //fastq_data = fastq_path.join(biosino_download.out.fastq_files, 
    //  remainder: true)
    // Trim adapters and low quality reads
    bbduk(fastq_path.concat(getFastq.out.public_fastq))
    // FastQC on trimmed reads
    fastqQCTrim(bbduk.out.trimmed_qc, "trimmed")
    // Align sequences to host genome, extract unalgined reads, align these
    // to the EBV genome, remove PCR duplicates, and extract the EBV reads
    hostRemoval(bbduk.out.trimmed_reads, human_index, ebv_index)
    // Split SAM file to get Host and Organism BAM files
    splitSam(hostRemoval.out.sam_path)
    // Extract Host reads
    extractHostReads(splitSam.out.host_bam_path)
    // Extract Organism reads
    extractOrgReads(splitSam.out.org_bam_path, ebv_index)
    // FastQC on EBV specific reads
    fastQCEbv(extractOrgReads.out.org_fastqc_input, "host_removed")
    // Repair the deduplicated EBV reads to appear as ordered pairs
    repairFastq(extractOrgReads.out.org_fastq_files)
    // Run spades on the EBV reads
    spades(repairFastq.out.repaired_fastq)
    // Run unicycler on the EBV reads
    unicycler(repairFastq.out.repaired_fastq)
    abacasSpades(spades.out.contigs.filter{ it.size()>0 }, "spades", ebv_index)
    abacasUnicycler(unicycler.out.contigs.filter{ it.size()>0 }, "unicycler", ebv_index)
    prokkaSpades(abacasSpades.out.abacas_contigs, "spades")
    prokkaUnicycler(abacasUnicycler.out.abacas_contigs, "unicycler")
    // Run quast on the EBV assembly
    quastSpades(abacasSpades.out.abacas_contigs.toSortedList(), ebv_masked_index, "spades")
    quastUnicycler(abacasUnicycler.out.abacas_contigs.toSortedList(), ebv_masked_index, "unicycler")
    haplotypeCaller(extractOrgReads.out.org_bam_path, ebv_masked_index)
    cnnScoreVariants(haplotypeCaller.out.vcf_paths, ebv_masked_index)
    // Annotate VCFs from simplex and duplex datasets
    annotateVCFs(cnnScoreVariants.out.annotated_vcf_paths.concat(duplex_vcf))
    // Align WGS reads to EBV genome and call variations
    bbdukWGS(wgs_path)
    alignEBV(bbdukWGS.out.trimmed_reads, ebv_index)
    haplotypeCallerWGS(alignEBV.out.bam_path, ebv_masked_index)
    cnnScoreVariantsWGS(haplotypeCallerWGS.out.vcf_paths, ebv_masked_index)
    annotateVCFsWGS(cnnScoreVariantsWGS.out.annotated_vcf_paths)
    createAlternateReferenceWGS(annotateVCFsWGS.out.snpeff_vcf_paths, ebv_masked_reference, ebv_findex, ebv_dict)
    // Perform gene annotation on all public, simplex and duplex datasets
    createAlternateReference(duplex_list, duplex_vcf.toSortedList(), duplex_tabix.toSortedList(), ebv_masked_reference, ebv_findex, ebv_dict)
    repeat_masker(abacasUnicycler.out.abacas_contigs.concat(createAlternateReference.out.alt_ref, public_fasta, createAlternateReferenceWGS.out.alt_ref), repeat_table)
    prokkaCombined(abacasUnicycler.out.abacas_contigs.concat(createAlternateReference.out.alt_ref, public_fasta, createAlternateReferenceWGS.out.alt_ref), "combined")
    find_epitope(repeat_masker.out.masked_fasta, epitope_db)
    quastCombined(repeat_masker.out.masked_fasta.toSortedList(), ebv_masked_index, "combined")
    selectGenomes(repeat_masker.out.masked_fasta.toSortedList(), quastCombined.out.report)
    variant_call_muscle(repeat_masker.out.masked_fasta, ebv_masked_reference)
    parsnp(selectGenomes.out.filtered_genomes, ebv_masked_reference)

}



workflow ebvAssemblyLite {
  ed_data = Channel.fromFilePairs(params.ed_data, size: 4)
  elshafa_data = Channel.fromFilePairs(params.elshafa_data, size: 2)
  fastq_path = ed_data.concat(elshafa_data)
  public_list = Channel.from(params.public_list)
  adapter_path = Channel.fromPath(params.adapters)
  public_fasta = Channel.fromPath(params.assembled_genome )
  out_path = Channel.fromPath(params.out_path)
  ebv_reference = Channel.fromPath(params.references.organism.ebv_1_ref + '/genome.fa')
  ebv_masked_reference = Channel.fromPath(params.references.masked.ebv_1_ref + '/genome.fa')
  human_reference = Channel.fromPath(params.references.host.sans_ebv + '/genome.fa')
  ebv_index = Channel.fromPath(params.references.organism.ebv_1_ref)
  ebv_masked_index = Channel.fromPath(params.references.masked.ebv_1_ref)
  human_index = Channel.fromPath(params.references.host.sans_ebv)
  multi_config = Channel.fromPath(params.multiqc_config)
  duplex_vcf = Channel.fromPath(params.duplex_vcfs)
  duplex_tabix = Channel.fromPath(params.duplex_tabix)
  duplex_list = Channel.from(params.duplex_list)
  ebv_findex = Channel.fromPath(params.references.masked.ebv_1_ref + '/genome.fa.fai')
  ebv_dict = Channel.fromPath(params.references.masked.ebv_1_ref + '/genome.dict')
  out_path = Channel.fromPath(params.out_path)
  repeat_table = Channel.fromPath(params.ebv_repeats)
  wgs_path = Channel.fromFilePairs(params.biosino_data,  size: 2)
  epitope_db = Channel.fromPath(params.ebv_epitope_db)
  aligned_ebv = Channel.fromFilePairs(params.biosino_alignments)
  study_metadata = Channel.fromPath(params.study_metadata)
  main:
    //getFastq(public_list)
    //Run FastQC on raw read
    fastqQCRaw(fastq_path, "raw")
    // Download biosino data
    //biosino_download(biosino_data.splitCsv(header: true, sep: '\t'))
    //fastqQCPublic(biosino_download.out.fastq_files.groupTuple(), "public")
    // Merge fastq channels
    //fastq_data = fastq_path.join(biosino_download.out.fastq_files, 
    //  remainder: true)
    // Trim adapters and low quality reads
    bbduk(fastq_path)
    // FastQC on trimmed reads
    fastqQCTrim(bbduk.out.trimmed_qc, "trimmed")
    // Align sequences to host genome, extract unalgined reads, align these
    // to the EBV genome, remove PCR duplicates, and extract the EBV reads
    hostRemoval(bbduk.out.trimmed_reads, human_index, ebv_index)
    // Split SAM file to get Host and Organism BAM files
    splitSam(hostRemoval.out.sam_path)
    // Extract Host reads
    extractHostReads(splitSam.out.host_bam_path)
    // Extract Organism reads
    extractOrgReads(splitSam.out.org_bam_path, ebv_index)
    // FastQC on EBV specific reads
    fastQCEbv(extractOrgReads.out.org_fastqc_input, "host_removed")
    // Repair the deduplicated EBV reads to appear as ordered pairs
    repairFastq(extractOrgReads.out.org_fastq_files)
    // Run spades on the EBV reads
    spades(repairFastq.out.repaired_fastq)
    // Run unicycler on the EBV reads
    unicycler(repairFastq.out.repaired_fastq)
    abacasSpades(spades.out.contigs.filter{ it.size()>0 }, "spades", ebv_index)
    abacasUnicycler(unicycler.out.contigs.filter{ it.size()>0 }, "unicycler", ebv_index)
    prokkaSpades(abacasSpades.out.abacas_contigs, "spades")
    prokkaUnicycler(abacasUnicycler.out.abacas_contigs, "unicycler")
    // Run quast on the EBV assembly
    quastSpades(abacasSpades.out.abacas_contigs.toSortedList(), ebv_masked_index, "spades")
    quastUnicycler(abacasUnicycler.out.abacas_contigs.toSortedList(), ebv_masked_index, "unicycler")
    haplotypeCaller(extractOrgReads.out.org_bam_path, ebv_masked_index)
    cnnScoreVariants(haplotypeCaller.out.vcf_paths, ebv_masked_index)
    // Annotate VCFs from simplex and duplex datasets
    annotateVCFs(cnnScoreVariants.out.annotated_vcf_paths.concat(duplex_vcf))
    // Align WGS reads to EBV genome and call variations
    haplotypeCallerWGS(aligned_ebv.combine(ebv_masked_index))
    cnnScoreVariantsWGS(haplotypeCallerWGS.out.vcf_paths, ebv_masked_index)
    annotateVCFsWGS(cnnScoreVariantsWGS.out.annotated_vcf_paths)
    createAlternateReferenceWGS(annotateVCFsWGS.out.snpeff_vcf_paths, ebv_masked_reference, ebv_findex, ebv_dict)
    // Perform gene annotation on all public, simplex and duplex datasets
    createAlternateReference(duplex_list, duplex_vcf.toSortedList(), duplex_tabix.toSortedList(), ebv_masked_reference, ebv_findex, ebv_dict)
    repeat_masker(abacasUnicycler.out.abacas_contigs.concat(createAlternateReference.out.alt_ref, public_fasta, createAlternateReferenceWGS.out.alt_ref), repeat_table)
    quastCombined(repeat_masker.out.masked_fasta.toSortedList(), ebv_masked_index, "combined")
    selectGenomes(repeat_masker.out.masked_fasta.toSortedList(), quastCombined.out.report, study_metadata)
    prokkaCombined(selectGenomes.out.fasta_files.collect(), "combined")
    find_epitope(selectGenomes.out.fasta_files.collect(), epitope_db)
    variant_call_muscle(selectGenomes.out.fasta_files.collect(), ebv_masked_reference)
    parsnp(selectGenomes.out.filtered_genomes, ebv_masked_reference)

}

workflow prepNetMHCPan {
  epitope_db = Channel.fromPath(params.ebv_epitope_db)
  fasta_paths = Channel.fromPath(params.masked_genomes)
  out_path = Channel.fromPath(params.out_path)
  main:
    find_epitope(fasta_paths, epitope_db)
} 

workflow callAllVariants {
  fasta_paths = Channel.fromPath(params.masked_genomes)
  ebv_masked_reference = Channel.fromPath(params.references.masked.ebv_1_ref + '/genome.fa')
  main:
    variant_call_muscle(fasta_paths, ebv_masked_reference)
}

workflow createConsensusFasta {
  vcf_path = Channel.fromPath(params.twinstrand_vcf)
  tabix_path = Channel.fromPath(params.twinstrand_tabix)
  sample_list = Channel.from(params.sample_list)
  reference_path = Channel.fromPath(params.ebv_1_ref + '/genome.fa')
  reference_index = Channel.fromPath(params.ebv_1_ref + '/genome.fa.fai')
  reference_dict = Channel.fromPath(params.ebv_1_ref + '/genome.dict')
  out_path = Channel.fromPath(params.out_path)
  main:
    create_alternate_reference(sample_list, vcf_path.toSortedList(), tabix_path.toSortedList(), reference_path, reference_index, reference_dict)

}

workflow annotateAllVCFs {
  vcf_paths = Channel.fromPath(params.vcf_paths)
  main:
    annotateVCFs(vcf_paths)
}

workflow alignAndCoverage {
  fastq_path = Channel.fromFilePairs(params.fastq_path, size: 2)
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
    //generateFastq(simulated_reads, org_reference_path)
    bbduk(fastq_path.combine(adapter_path))
    bowtie_Align(bbduk.out.trimmed_reads.combine(index_path))
    sam_to_bam(bowtie_Align.out.sam_path)
    bam_sort(sam_to_bam.out.bam_files)
    keep_unaligned(bam_sort.out.sorted_bam_file)
    align_ebv(bbduk.out.trimmed_reads.combine(org_index_path))
    align_filtered_ebv(keep_unaligned.out.filtered_reads.combine(org_index_path))
    haplotype_caller(align_filtered_ebv.out.sam_path.combine(org_reference_path))
    cnn_score_variants(haplotype_caller.out.vcf_paths.combine(org_reference_path))
    annotate_vcfs(cnn_score_variants.out.annotated_vcf_paths)
    bbnorm(keep_unaligned.out.filtered_reads)
    align_sampled_ebv(bbnorm.out.sampled_reads.combine(org_index_path))
    spades(bbnorm.out.sampled_reads.combine(org_reference_path))
    quast("Spades", spades.out.contigs.toSortedList(), org_reference_path)
    unicycler(bbnorm.out.sampled_reads)
    quast_uniclyer("Unicycler", unicycler.out.contigs.toSortedList().unique(), org_reference_path)
    abacas(unicycler.out.unicycler.combine(org_reference_path))
    quast_abacas("Abacas", abacas.out.contigs.toSortedList(), org_reference_path)    
    repeat_masker(abacas.out.contigs.concat(public_contigs), ebv_repeat)
    quast_masker("RepeatMasker", repeat_masker.out.contigs.toSortedList(), org_reference_path)
    mpileup(align_sampled_ebv.out.sam_path.combine(org_reference_path))
    prokka(unicycler.out.contigs,"Unicycler")
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
  public_fastq = Channel.fromSRA(params.public_list, protocol: 'http', apiKey: params.api_key)
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

workflow annotateAllGenomes {
  fasta_path = Channel.fromPath(params.fasta_path)
  main:
    prokka(fasta_path,"annotated")

}