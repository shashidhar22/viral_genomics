# Viral genomics

This pipeline is primarily used to call variant and assemble short-read Illumina 
sequencing data from viral samples. The pipeline is developed using Nextflow. 
The pipelines follow the Nextflow DSL2 guidelines, and each process has a 
singularity container associated with it. Except for `genome_annotation.nf` 
which uses a custom conda environment that is included with the repo. 

## ebvAssemble

The primary workflow is `ebvAssembly`. The inputs for the workflow are read from
`main.yaml`. We use the `hg38` human genome with the `EBV` contig removed to 
remove host human sequences from our sequencing samples. The `EBV` reference 
`NC_007605.1 Human gammaherpesvirus 4, complete genome`. The workflow takes the 
following steps:

1. QC on raw reads (`fastqQCRaw`)
2. Trim adapters and low quality reads (`bbduk`)
3. QC on trimmed reads (`fastqQCTrim`)
4. Remove host reads and extract EBV reads (`hostReads`)
5. QC on EBV specific reads (`fastQCEbv`)
6. Repair read pair orders in EBV FASTQ files (`repairFastq`)
7. Run SPAdes assembly (`spades`)
8. Run Unicycler assembly (`unicycler`)
9. Run ABACAS on both SPAdes and Unicycler assemblies (`abacasSpades`, `abacasUnicyler`)
10. Identify genes from assemblies (`prokkaSpades`, `prokkaUnicycler`)
11. Generate MultiQC reports (`runMultiQC`)
12. Generate Quast reports (`quastSpades`, `quastUnicycler`)
13. Call variants and annotate variants using GATK and SnpEff (`haplotypeCaller`, `cnnScoreVariants`, `annotateVCFs`)

A run script `run.sh` is included in the repertoire to execute the workflow. The
user would need to edit the following variables in `run.sh` before running the 
script
1. `NXF_CONFIG`
2. `WORKFLOW_REPO`
3. `QUEUE`
4. `path_to_work_dir`
5. `path_to_yaml.yaml`

The provided `main.yaml` serves as a template. Reference genomes used were 
obtained here:
1. Human genome: http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/NCBI/GRCh38/Homo_sapiens_NCBI_GRCh38.tar.gz
2. EBV type 1 reference genome: https://www.ncbi.nlm.nih.gov/nuccore/NC_007605.1/
3. EBV type 2 reference genome: https://www.ncbi.nlm.nih.gov/nuccore/NC_009334.1/
4. GATK databases: https://github.com/bahlolab/bioinfotools/blob/master/GATK/resource_bundle.md


Note: We generated the reference `sans_ebv` mentioned in `main.yaml` by removing
the `EBV` contig from the `hg38` genome file

## Setting up SnpEff databases:

1. Change directory to `viral_genomics` folder
2. Install the `ebv_enktl` environment
   
   ```{sh}
   conda env create -f snpeff.yml
   ```

3. Change the path to `data.dir` in `snpEff.config` to the `ebv` folder in `envs`   
4. Copy the `snpEff.config` file to the correct path
   
   ```{sh}
   cp envs/snpEff.config ~/.conda/envs/ebv_enktl/share/snpeff-5.0-1/
   ```
