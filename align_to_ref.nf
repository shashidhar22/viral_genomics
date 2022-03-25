// Container versions
container__bwa = "quay.io/fhcrc-microbiome/bwa:bwa.0.7.17__bcw.0.3.0I"

// Optional parameters
// Minimum alignment score which will result in a set of reads being retained
// after alignment to the reference genome
params.min_align_score = 20


// PROCESSES

// Process to build a BWA index for a genome
process bwa_index {
    container "${container__bwa}"

    input:
        file fasta

    output:
        file "reference.tar.gz"

"""#!/bin/bash
set -e
# Build the BWA index
bwa index ${fasta}
# Make the tarball
tar cvf reference.tar "${fasta}*"
# Compress it
gzip reference.tar
"""
}

// Process to retain reads which align to a reference sequence
process bwa_filter {
    container "${container__bwa}"

    input:
        tuple file(index_tgz), val(sample_name), file(R1), file(R2)

    output:
        tuple val(sample_name), file("${R1}.filtered.fq.gz"), file("${R2}.filtered.fq.gz")


"""#!/bin/bash
# Stop execution if any of the individual commands have a non-zero exit status
set -e
# Get the name of the indexed genome
bwa_index_fn=\$(tar -ztvf ${index_tgz} | head -1 | sed \'s/.* //\')
bwa_index_prefix=\${bwa_index_fn%.*}
echo BWA index prefix is \${bwa_index_prefix}
# Extract the tarball containing the indexed genome
echo Extracting BWA index
mkdir -p bwa_index/ 
tar -I pigz -xf ${index_tgz} -C bwa_index/
echo Files in index directory: 
ls -l -h bwa_index 
echo Running BWA 
bwa mem -t ${task.cpus} \
	-T ${params.min_align_score} \
	-o alignment.sam \
	bwa_index/\$bwa_index_prefix \
	${R1} ${R2}
# The line below will raise an error if the file is empty
echo Checking if alignment is empty
[[ -s alignment.sam ]]
# In the comand below, the "-F 12" flag will instruct samtools to filter
# out any pairs in which both reads are unmapped. This will retain any
# pairs in which either read is mapped.
# For an intuitive explanation of these numeric flags, see 
# https://broadinstitute.github.io/picard/explain-flags.html
echo Extracting Unaligned Pairs 
samtools fastq alignment.sam \
	--threads ${task.cpus} -F 12 \
	-1 ${R1}.filtered.fq.gz -2 ${R2}.filtered.fq.gz
echo Done
"""
}


// WORKFLOW

workflow {

    reference_fasta = Channel.fromPath(params.ref_fasta)
    input_reads = Channel.fromFilePairs(params.fastq_paths, size: 4)

    main:
        bwa_index(reference_fasta)
        bwa_filter(bwa_index.out.combine(input_reads)

}