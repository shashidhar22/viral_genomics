nextflow.enable.dsl=2

// Call variants from reads aligned to a reference genome using GATK 
// HaplotypeCaller. Prior to variant calling read group information is added
// to the reads using the GATK AddOrReplaceReadGroups tool.
//
// PROCESS INPUTS:
// - Tuple with sample name, path sorted bam file, the bam index file, and the
//  reference genome fasta file.
//  (Source: align_filtered_ebv.out.sam_path.combine(org_reference_path))
//
// PROCESS OUTPUTS:
// - Raw VCF files for each sample with tabix index file
// EMITTED CHANNELS:
// - vcf_paths :- Tuple with sample name, path to raw vcf file, and the tabix
// file
//
// NOTE: 
//
// TODO: 
process haplotype_caller {
  container "$params.gatk"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/variant_calls/${sample}/", mode : "copy"
  input:
    each path(bam_path)
    path reference_path
  output:
    path "${sample}.vcf.gz*" , emit: vcf_paths
    
  script:
    def bam_file = bam_path[1]
    sample =  bam_file.getSimpleName().replaceAll(/_ebv_dedup/, "")
    """
    gatk AddOrReplaceReadGroups \
      -I ${bam_file} \
      -O ${sample}_rg.bam \
      -LB ebv \
      -PL ILLUMINA \
      -PU NA \
      -SM ${sample} \
      --CREATE_INDEX  true
    gatk --java-options "-Xmx4g" HaplotypeCaller \
      -R ${reference_path}/genome.fa \
      -I ${sample}_rg.bam \
      -ploidy 1 \
      --output ${sample}.vcf.gz  > stdout.log 2> stderr.log 
    """
}


process haplotype_caller_wgs {
  container "$params.gatk"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/variant_calls/${sample}/", mode : "copy"
  input:
    tuple val(sample), path(bam_path), path(reference)
  output:
    path "${sample}.vcf.gz*" , emit: vcf_paths
  script:
    def bam_file = bam_path[1]
    """
    gatk AddOrReplaceReadGroups \
      -I ${bam_file} \
      -O ${sample}_rg.bam \
      -LB ebv \
      -PL ILLUMINA \
      -PU NA \
      -SM ${sample} \
      --CREATE_INDEX  true
    gatk --java-options "-Xmx4g" HaplotypeCaller \
      -R ${reference}/genome.fa \
      -I ${sample}_rg.bam \
      -ploidy 1 \
      --output ${sample}.vcf.gz  > stdout.log 2> stderr.log 
    """
}

// Score variants using a convolution neural network model trained on the 
// variant calling data from the previous step and the reference context, 
// providing a score for each variant. Here we use the CNNScoreVariants toolkit
// from the GATK to score the variants.
//
// PROCESS INPUTS:
// - Tuple with sample name, path to VCF file, the tabix index file, and the
//  reference genome fasta file.
//  (Source: haplotype_caller.out.vcf_paths.combine(org_reference_path))
//
// PROCESS OUTPUTS:
// - Variants scored using the CNNScoreVariants toolkit.
// EMITTED CHANNELS:
// - annotated_vcf_paths :- Tuple with sample name, path to scored vcf file
//
// NOTE: 
//
// TODO: 
process cnn_score_variants {
  container "$params.gatk"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.out_path/variant_calls/${sample}/", mode : "copy"
  input:
    each path(vcf_path)
    path reference_path
  output:
    path "${sample}_annotated.vcf", emit: annotated_vcf_paths
  script:
    def vcf_file = vcf_path[0]
    sample = vcf_file.getSimpleName()
    """
    mkdir temp_dir
    gatk CNNScoreVariants \
      -V ${vcf_file} \
      -R ${reference_path}/genome.fa \
      --tmp-dir temp_dir \
      -O ${sample}_annotated.vcf > stdout.log 2> stderr.log
    """
}

// Annotate variant using SnpEff. For this situation, a ebv specific DB was 
// created using the SnpEff command line tool. To run the module, you will 
// need to have the conda environment ebv_env.yml installed.
//
// PROCESS INPUTS:
// - Tuple with sample name, path to VCF file,
//  (Source: cnn_score_variants.out.annotated_vcf_paths)
//
// PROCESS OUTPUTS:
// - Variants annotated with functional information using SnpEff, along with
// a HTML summary of the annotations.
// EMITTED CHANNELS:
// - snpeff_vcf_paths :- Tuple with sample name, path to annotated vcf file,
// and the html summary file
//
// NOTE: 
//
// TODO: Create a conda environment for SnpEff (ebv_env.yml)
process annotate_vcfs {
  conda "$params.conda.enktl"
  errorStrategy 'retry'
  time '1d 6h'
  label 'local'
  cache true
  publishDir "$params.out_path/variant_calls/${sample}/", mode : "copy"
  input:  
    each path(vcf_path) 
  output:
    path "${sample}_snpeff.vcf" , emit: snpeff_vcf_paths
    path "${sample}_summary.html" , emit: snpeff_vcf_summary
    path "${sample}_snpeff.tsv", emit: snpeff_vcf_tsv
  script:
    sample = vcf_path.getSimpleName().replaceAll("_annotated", "")
    """
    snpEff -v ebv -s ${sample}_summary.html ${vcf_path} > ${sample}_snpeff.vcf
    cat ${sample}_snpeff.vcf | ~/.conda/pkgs/snpsift-4.3.1t-hdfd78af_3/share/snpsift-4.3.1t-3/scripts/vcfEffOnePerLine.pl  > ${sample}_oneperline.vcf 
    SnpSift extractFields ${sample}_oneperline.vcf CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].HGVS_P" > ${sample}_snpeff.tsv 
    """
}

// Call variants from the host genome using GATK HaplotypeCaller.
//
// PROCESS INPUTS:
// - Tuple with sample name, path to sorted BAM file, the BAM index file, and
//  the host reference genome 
//  (Source: bam_sort.out.sorted_bam_file.combine(varcall_host_path))
//
// PROCESS OUTPUTS:
// - GVCF files for each sample with tabix index file
// EMITTED CHANNELS:
// - gvcf_paths :- Path to GVCF file
// - gvcf_index_paths :- Path to tabix index file
//
// NOTE: 
//
// TODO: 
process call_host_variants {
  container "$params.gatk"
  errorStrategy 'retry'
  time '1d 6h'
  label 'mid_mem'
  publishDir "$params.out_path/host_variant_calls/", mode : "copy"
  input:
    each path(bam_file)
    path reference_path
  output:
    path "${sample}.g.vcf.gz", emit: gvcf_paths
    path "${sample}.g.vcf.gz.tbi", emit: gvcf_index_paths
  script:
    def bam_path = bam_file[0]
    def bam_index_path = bam_file[1]
    def sample =  ${bam_path}.replaceAll(/_ebv_dedup.bam/, "")
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller \
      -R ${reference_path}/genome.fa \
      -I ${bam_file} \
      -O ${sample}.g.vcf.gz \
      -ERC GVCF
    """
}

// Generate a config file mapping sample names to gVCF files and prepare the 
// input for the consolidate gVCF step.
//
// PROCESS INPUTS:
// - Path to the gVCF file 
//  (Source: call_host_variants.out.gvcf_paths)
// - Path to the tabix index file
//  (Source: call_host_variants.out.gvcf_index_paths)
//
// PROCESS OUTPUTS:
// - Path to file mapping sample names to gVCF files
// EMITTED CHANNELS:
// - gvcf_map_file :- Path to map file
//
// NOTE: 
//
// TODO: 
process generate_gvcf_table {
  conda "$params.conda.enktl"
  errorStrategy 'retry'
  time '1d 6h'
  label 'mid_mem'

  publishDir "$params.out_path/variant_calls/", mode : "copy"
  input:
    path gvcf_paths
    path gvcf_index_paths
  output:
    path "gvcf_map_file.tsv", emit: gvcf_map_file
  script:
    """
    Rscript $moduleDir/rscripts/map_gvcfs.R -p ./ -o gvcf_map_file.tsv 
    """
}

// Generate a config file mapping sample names to gVCF files and prepare the 
// input for the consolidate gVCF step.
//
// PROCESS INPUTS:
// - Path to the gVCF file 
//  (Source: call_host_variants.out.gvcf_paths)
// - Path to the tabix index file
//  (Source: call_host_variants.out.gvcf_index_paths)
//
// PROCESS OUTPUTS:
// - Path to file mapping sample names to gVCF files
// EMITTED CHANNELS:
// - gvcf_map_file :- Path to map file
//
// NOTE: 
//
// TODO: 
process consolidate_gvcfs {
  conda "$params.conda.enktl"
  errorStrategy 'retry'
  time '1d 6h'
  label 'mid_mem'

  publishDir "$params.out_path/host_variant_calls/", mode : "copy"
  input:
    path gvcf_paths 
    path gvcf_index_paths
    path gvcf_map_file
  output:
    path "gvcf_database", emit: gvcf_database
  script:
    """
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      --genomicsdb-workspace-path gvcf_database --batch-size 0 \
      -L 20 --sample-name-map gvcf_map_file.tsv \
      --reader-threads 4 > stdout.log 2> stderr.log  
    """
}

process genotype_gvcfs {
  container "$params.gatk"
  errorStrategy 'retry'
  time '1d 6h'
  label 'mid_mem'
  publishDir "$params.out_path/host_variant_calls/", mode : "copy"
  input:
    path gvcf_database
    path reference_path
  output:
    path "host_variant_calls.vcf.gz", emit : joint_vcf
    path "host_variant_calls.vcf.gz.tbi", emit : joint_vcf_index
  script:
    """
    mkdir temp_dir
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
      -R ${reference_path}/genome.fa \
      -V gendb://gvcf_database \
      --tmp-dir temp_dir \
      -O host_variant_calls.vcf.gz > stdout.log 2> stderr.log  
    """
}

process variant_scoring {
  container "$params.gatk"
  errorStrategy 'retry'
  time '1d 6h'
  label 'mid_mem'
  publishDir "$params.out_path/host_variant_calls/", mode : "copy"  
  input:
    path joint_vcf
    path joint_vcf_index
    path reference_path
  output:
    path "host_variant_calls.vcf.gz", emit : annotated_vcf
    path "host_variant_calls.vcf.gz.tbi", emit : annotated_vcf_index
  script:  
    """ 
    mkdir temp_dir
    gatk VariantRecalibrator \
      --tmp-dir temp_dir \
      -R ${reference_path}/genome.fa \
      -V host_variant_calls_scored.vcf.gz \
      --resource hapmap,known=false,training=true,truth=true,prior=15.0:$params.variants.hapmap \
      --resource mills,known=false,training=true,truth=true,prior=12:$params.variants.mills \
      --resource axiomPoly,known=false,training=true,truth=false,prior=10:$params.variants.axiom \
      --resource omni,known=false,training=true,truth=false,prior=12.0:$params.variants.1000gomni \
      --resource 1000G,known=false,training=true,truth=false,prior=10.0:$params.variants.1000gphase3 \
      --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$params.variants.dbsnp \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
      -mode BOTH \
      -O host_variant_calls.recal \
      --tranches-file host_variant_calls.tranches 
    gatk ApplyVQSR \
      --tmp-dir temp_dir \
      -R ${reference_path}/genome.fa \
      -V host_variant_calls.vcf.gz \
      -O host_variant_calls_scored.vcf.gz \
      --ts_filter_level 99.0 \
      --tranches-file host_variant_calls.tranches \
      --recal-file host_variant_calls.recal \
      -mode BOTH
    """
}



process create_alternate_reference {
  container "$params.gatk"
  errorStrategy 'retry'
  time '1h'
  label 'mid_mem'
  publishDir "$params.out_path/alternate_reference/", mode : "copy"  
  input:
    each sample 
    path vcf_data
    path tabix_data
    path ref_path
    path ref_index
    path ref_dict

  output:
    path "${sample}_consensus.fasta", emit: alt_ref
  script:
    """
    gatk FastaAlternateReferenceMaker -R ${ref_path} -O ${sample}.fasta -V ${sample}.consensus.variant-calls.vcf.gz
    awk '/^>/ {printf(">%s\\n","${sample}");next;} {print}' ${sample}.fasta > ${sample}_consensus.fasta 
    """

}

process create_alternate_reference_wgs {
  container "$params.gatk"
  errorStrategy 'retry'
  time '1h'
  label 'mid_mem'
  publishDir "$params.out_path/alternate_reference/", mode : "copy"  
  input:
    each path(vcf_path)
    path ref_path
    path ref_index
    path ref_dict

  output:
    path "${sample}_consensus.fasta", emit: alt_ref
  script:
    sample = vcf_path.getSimpleName().replaceAll("_snpeff", "")
  
    """
    gatk IndexFeatureFile -I ${sample}_snpeff.vcf
    gatk FastaAlternateReferenceMaker -R ${ref_path} -O ${sample}.fasta -V ${sample}_snpeff.vcf
    awk '/^>/ {printf(">%s\\n","${sample}");next;} {print}' ${sample}.fasta > ${sample}_consensus.fasta 
    """

}




process variant_call_muscle {
  conda "$params.conda.enktl"
  errorStrategy 'retry'
  time '1h'
  label 'mid_mem'
  publishDir "$params.out_path/whole_genome_variant_calls/", mode : "copy", pattern : "*.tsv"
  publishDir "$params.out_path/whole_genome_variant_calls_vcfs/", mode : "copy", pattern : "*.vcf" 
  input:
    each path(fasta_file)
    path ref_path

  output:
    path "${sample}_snpeff.tsv", emit: alt_ref
    path "${sample}.vcf", emit: all_vcf

  script:
    sample = fasta_file.getSimpleName().replaceAll(".\\d?.fasta.masked|.\\d_consesus.fasta.masked|_S\\d+_abacas.fasta.masked|_2300_abacas.fasta.masked|_consensus/.fasta.masked", "")
    """
    nucmer -p ${sample} ${ref_path} ${fasta_file}
    show-snps -T ${sample}.delta > ${sample}.snps
    /home/sravisha/software/all2vcf/all2vcf mummer --no-Ns --snps ${sample}.snps --reference ${ref_path} --input-header --output-header > ${sample}.vcf 
    snpEff -Xmx4g -v ebv -s ${sample}_summary.html ${sample}.vcf > ${sample}_snpeff.vcf
    cat ${sample}_snpeff.vcf | ~/.conda/pkgs/snpsift-4.3.1t-hdfd78af_3/share/snpsift-4.3.1t-3/scripts/vcfEffOnePerLine.pl  > ${sample}_oneperline.vcf 
    SnpSift extractFields ${sample}_oneperline.vcf CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].HGVS_P" > ${sample}_snpeff.tsv 
    """

}

process variant_call_genome {
  container "$params.gatk"
  errorStrategy 'retry'
  time '1h'
  label 'mid_mem'
  publishDir "$params.out_path/variant_call_genome/", mode : "copy"  
  input:
    each path(vcf_path)
    path ref_path
    path ref_index
    path ref_dict

  output:
    path "${sample}_consensus.fasta", emit: alt_ref
  script:
    sample = vcf_path.getSimpleName().replaceAll("_sorted_snpeff", "")
  
    """
    gatk IndexFeatureFile -I ${sample}_snpeff.vcf
    gatk FastaAlternateReferenceMaker -R ${ref_path} -O ${sample}.fasta -V ${sample}_snpeff.vcf
    awk '/^>/ {printf(">%s\\n","${sample}");next;} {print}' ${sample}.fasta > ${sample}_consensus.fasta 
    """

}