# mutect-nf
## Mutect pipeline for somatic variant calling with Nextflow
[![CircleCI](https://circleci.com/gh/IARCbioinfo/RNAseq-nf/tree/master.svg?style=svg)](https://circleci.com/gh/nalcala/mutect-nf/tree/master)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/nalcala/mutect-nf/)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4333)

![workflow](mutectseqpipeline.png?raw=true "Scheme of calling Workflow")

## Description
Nextflow pipeline for somatic variant calling with mutect with Mutect1 or 2, gatk3 or gatk4

## Dependencies
1. Nextflow: for common installation procedures see the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository.
2. [Mutect](https://github.com/broadinstitute/mutect) and its dependencies (Java 1.7 and Maven 3.0+), or [gatk4](https://github.com/broadinstitute/gatk) that now includes Mutect2
3. [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html) and move the executable file in your path.
4. [python](https://www.python.org/) and package [pysam](https://github.com/pysam-developers/pysam)
5. [bedops](https://github.com/bedops/bedops)

**A conda receipe, and docker and singularity containers are available with all the tools needed to run the pipeline (see "Usage")**

### GATK4
With GATK4, a list of known_snps can be provided to mutect2 to improve the variant classification, for example file [af-only-gnomad.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz) from the bundle best practices from the broad institute [GATK somatic calling bundle](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38/).

### estimate contamination
When the estimate contamination mode is chosen, one needs to provide a list of known snps; we recommend the file [small_exac_common_3.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz) from the best practices broad institute bundle.

## Input 
 | Type      | Description     |
  |-----------|---------------|
  |--bam_folder    | a folder with tumor and normal bam files |
  |--tumor_bam_folder | a folder with tumor bam files |
  |--normal_bam_folder | a folder with normal bam files |
  |--tn_file |  input tabulation-separated values file with columns SM (sample name), RG (read group), pair1 (first fastq pair file), and pair2 (second fastq pair file) |
  
  Note that there are three input methods: single bam_folder, separated tumor_bam_folder and normal_bam_folder, and tn_file. The bam_folder method is the easiest and assumes that both normal and tumor bam files are in this folder, and it uses parameters suffix_tumor and suffix_normal to detect them (the rest of the file name needs to be identical); the separated tumor_bam_folder and normal_bam_folder method also uses the suffix_tumor and suffix_normal to match the samples. 
  
  The tn_file method uses a tabulation-separated values format file with columns sample, tumor, and normal (in any order); it does not use parameters suffix_tumor and suffix_normal and does not require file names to match. When the genotype mode is active, additional columns are expected: preproc, specifying if preprocessing of RNA-seq bam file is required (yes or no) and vcf, indicating the location of the vcf file containing the alleles to genotype. The tn_file method is necessary for joint multi-sample calling, in which case the sample name is used to group files.

## Parameters

* #### Mandatory
| Name | Example value | Description |
|-----------|--------------:|-------------| 
|--ref | ref.fa | reference genome fasta file |
|--bed   |  gene.bed | bed file with genes for RESeQC | 


* #### Optional

```mutect_path```, ```--dbsnp```, ```--cosmic```

| Name | Default value | Description |
|-----------|--------------|-------------| 
|--cpu          | 4 | number of CPUs |
|--mem         | 8 | memory for mapping|
|--suffix_tumor      | \_T | suffix for tumor file|
|--suffix_normal      | \_N | suffix for matched normal file|
|--output_folder   | mutect_results | output folder for aligned BAMs|
|--known_snp |  "" | VCF file with known variants and frequency (e.g., from gnomad) |
|--mutect_args = ""
|--nsplit = 1
|--region = null
|--bed = null
|--java = "java"
|--snp_contam = "NO_FILE"
|--cosmic = ""
|--mutect_jar = null
|--mutect2_jar = null
|--gatk_version= "4"
|--PON = null

* #### Flags

| Name  | Description |
|-----------|-------------| 
|--help | print usage and optional parameters |
|--cutadapt | enable adapter and quality reads trimming before alignment|
|--sjtrim   | enable reads trimming at splice junctions | 
|--hisat2   | use hisat2 instead of STAR for mapping | 
|--recalibration  | perform quality score recalibration (GATK)|

estimate_contamination = null
genotype = null
RNAseq_preproc = null


## Usage
Nextflow seamlessly integrates with GitHub hosted code repositories:

`nextflow run iarcbioinfo/mutect-nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --bed mybedfile.bed --ref ref.fasta --mutect_jar mutect.jar`

If you specify option `--mutect2_jar` (GATK executable jar, which integrate mutect2) instead of `--mutect_jar`, the pipeline will automatically switched to mutect version 2.

#### Help section
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/mutect-nf --help
```
This shows details about optional and mandatory parameters provided by the user.  

#### BAM file format
The tumor bam file format must be (`sample` `suffix_tumor` `.bam`) with `suffix_tumor` as `_T` by default and customizable in input (`--suffix_tumor`). (e.g. `sample1_T.bam`)
The normal bam file format must be (`sample` `suffix_normal` `.bam`) with `suffix_normal` as `_N` by default and customizable in input (`--suffix_normal`). (e.g. `sample1_N.bam`).
BAI indexes have to be present in the same location than their BAM mates, with the extension `bam.bai`.

## Output 
  | Type      | Description     |
  |-----------|---------------|
  | BAM/file.bam    | BAM files of alignments or realignments |
  | BAM/file.bam.bai    | BAI files of alignments or realignments |
  | BAM/STAR.file.Chimeric.SJ.out.junction | STAR chimeric junction output |
  | BAM/STAR.file.SJ.out.tab | STAR junction tab output |
  | counts/file_count.txt                   | htseq-count output file  |
  | QC/multiqc_pretrim_report.html  | multiqc report before trimming | 
  | QC/multiqc_pretrim_report_data            | folder with data used to compute multiqc report before trimming |
  | QC/multiqc_posttrim_report.html      |     multiqc report before trimming | 
  | QC/multiqc_posttrim_report_data      |  folder with data used to compute multiqc report before trimming |
  | QC/adapter_trimming/file_{12}.fq.gz_trimming_report.txt | trim_galore report | 
  | QC/adapter_trimming/file_{12}_val_{12}_fastqc.zip | FastQC report after trimming | 
  | QC/alignment/STAR.file.Log.final.out, STAR.file.Log.out, STAR.file.Log.progress.out | STAR logs |
  | QC/bam/file_readdist.txt, file_clipping_profile\*, file_jun_stauration\*| RSeQC reports |
  | QC/fastq/file_{12}_pretrim_fastqc.zip | FastQC report before trimming | 
  | QC/BAM/BQSR/file_recal.table | table of scores before recalibration   |
  | QC/BAM/BQSR/file_post_recal.table   | table of scores after recalibration |
  | QC/BAM/BQSR/file_recalibration_plots.pdf   |  before/after recalibration plots   |
          
The output_folder directory contains three subfolders: BAM, counts, and QC


## Directed Acyclic Graph
[![DAG](dag.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/mutect-nf/blob/dev/dag.html)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------| 
  | Nicolas Alcala*    | AlcalaN@iarc.fr    | Developer to contact for support |
  | Tiffany Delhomme |  | Developer |
