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

#### Help section
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run mutect.nf --help
```
This shows details about optional and mandatory parameters provided by the user.  

#### BAM file format
The tumor bam file format must be (`sample` `suffix_tumor` `.bam`) with `suffix_tumor` as `_T` by default and customizable in input (`--suffix_tumor`). (e.g. `sample1_T.bam`)
The normal bam file format must be (`sample` `suffix_normal` `.bam`) with `suffix_normal` as `_N` by default and customizable in input (`--suffix_normal`). (e.g. `sample1_N.bam`).
BAI indexes have to be present in the same location than their BAM mates.

#### Global parameters
```mutect_path```, ```--dbsnp```, ```--cosmic```, and ```--ref``` are mandatory parameters but can be defined in your nextflow config file (```~/.nextflow/config``` or ```config``` in the working directory) and so not set as inputs.

The following is an example of config part defining this:
```bash
profiles {

        standard {
                params {
                   ref = '~/Documents/Data/references/hg19.fasta'
                   mutect = '~/Documents/Softs/MuTect/mutect-src/mutect/target/mutect-1.1.7.jar'
                   dbsnp = '~/Documents/Data/references/dbsnp_138.hg19_noMT.vcf'
                   cosmic = '~/Documents/Data/references/Cosmic_v73.hg19_noMT.vcf'
                }
        }
```
