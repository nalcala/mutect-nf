## mutect-nf
### Mutect pipeline on tumor-matched normal bam folder with Nextflow

#### Dependencies
1. Install [Mutect](https://github.com/broadinstitute/mutect) and its dependencies (Java 1.7 and Maven 3.0+).

2. Install [nextflow](http://www.nextflow.io/).

	```bash
	curl -fsSL get.nextflow.io | bash
	```
	And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
	```bash
	sudo mv nextflow /usr/local/bin
	```


#### Execution 
Either Nextflow is in your path and you have relevant rights on the script:   
`./mutect.nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --bed mybedfile.bed --ref ref.fasta`  
Or :  
`nextflow run mutect.nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --bed mybedfile.bed --ref ref.fasta`

#### Help
You can print the help manual by providing `--help` input in the execution command line.  
This shows details about optional and mandatory parameters provided by the user.  

#### BAM file format
The tumor bam file format must be (`sample` `suffix_tumor` `.bam`) with `suffix_tumor` as `_T` by default and customizable in input (`--suffix_tumor`). (e.g. `sample1_T.bam`)
The normal bam file format must be (`sample` `suffix_normal` `.bam`) with `suffix_normal` as `_N` by default and customizable in input (`--suffix_normal`). (e.g. `sample1_N.bam`).
BAI indexes have to be present in the same location than their BAM mates.
