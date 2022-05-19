nextflow.enable.dsl=2

// Trim reads using BBDuk to remove adapter sequences and lower quality bases
//
// PROCESS INPUTS:
// - Tuple with sample name, path to raw FASTQ file, and path to adapter file
// (Source: fastq_path.concat(generateFastq.out.sim_fastq).combine(adapter_path))
//
// PROCESS OUTPUTS:
// - Trimmed FASTQ files with adapter and low quality bases removed
// EMITTED CHANNELS:
// - trimmed_reads :- Channel containing sample name, and path to trimmed FASTQ 
// file
//
// NOTE: 
//
// TODO: 
process bbduk {
  errorStrategy 'retry'
  container "$params.bbmap"
  label 'high_mem'
  //publishDir "$params.out_path/trimmed_reads/", mode : "copy"
  input:
    tuple val(sample), path(fastq_path)
  output:
    path "*clean.fq", emit: trimmed_reads
    tuple val(sample), path("*clean.fq"), emit: trimmed_qc
    path "${sample}_bbduk_stats.txt", emit: trimmed_stats
  script:
  if ( fastq_path.size() == 2 )
    """
    bbduk.sh -Xmx400g in1=${fastq_path[0]} in2=${fastq_path[1]} \
    out1=${sample}_R1.clean.fq out2=${sample}_R2.clean.fq \
    ref=adapters qtrim=rl minlength=100 stats=${sample}_bbduk_stats.txt \
    ktrim=r k=23 mink=11 hdist=1 tpe=f tbo=t > stdout.log 2> stderr.log
    """
  else
    """
    zcat ${fastq_path[0]} ${fastq_path[2]} > ${sample}_R1.fastq 
    zcat ${fastq_path[1]} ${fastq_path[3]} > ${sample}_R2.fastq 
    bbduk.sh -Xmx400g in1=${sample}_R1.fastq in2=${sample}_R2.fastq \
    out1=${sample}_R1.clean.fq out2=${sample}_R2.clean.fq stats=${sample}_bbduk_stats.txt\
    ref=adapters qtrim=rl minlength=100 \
    ktrim=r k=23 mink=11 hdist=1 tpe=f tbo=t > stdout.log 2> stderr.log
    """
} 

// Normalize read depth across the genome using BBNorm to ensure all regions
// are equally sampled
//
// PROCESS INPUTS:
// - Tuple with sample name, path to FASTQ files containing organism specific 
//  reads (Source: keep_unaligned.out.filtered_reads)
//
// PROCESS OUTPUTS:
// - Depth normalized FASTQ files 
// EMITTED CHANNELS:
// - sampled_reads :- Channel containing sample name, and path to depth 
//  normalized FASTQ reads
//
// NOTE: 
//
// TODO: 
process bbnorm {
  errorStrategy 'retry'
  container "$params.bbmap"
  label 'mid_mem'
  publishDir "$params.out_path/sampled_reads/", mode : "copy"
  input:
    tuple val(sample), path(fastq_path)
  output:
    tuple val(sample), path("*_sampled_R*.fq"), emit: sampled_reads
  script:
  if ( fastq_path.size() == 2 )
    """
    mkdir -p temp_storage
    export TMPDIR=$workDir
    bbnorm.sh -Xmx20g in=${fastq_path[0]} in2=${fastq_path[1]} \
    out=${sample}_sampled_R1.fq out2=${sample}_sampled_R2.fq \
    target=1000 tmpdir=temp_dir 1> stdout.log 2> stderr.log
    """
  else
    """
    mkdir -p temp_storage
    export TMPDIR=$workDir
    zcat ${fastq_path[0]} ${fastq_path[2]} > ${sample}_R1.fastq 
    zcat ${fastq_path[1]} ${fastq_path[3]} > ${sample}_R2.fastq 
    bbnorm.sh -Xmx20g in=${sample}_R1.fastq in2=${sample}_R2.fastq \
    out=${sample}_sampled_R1.fq out2=${sample}_sampled_R2.fq \
    target=1000  tmpdir=temp_dir 1> stdout.log 2> stderr.log
    """
} 

