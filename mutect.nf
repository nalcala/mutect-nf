#! /usr/bin/env nextflow

// usage : ./mutect.nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --bed mybedfile.bed --ref ref.fasta --mutect_args " --force_output --force_alleles "

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '                 NEXTFLOW MUTECT                  '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run mutect.nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --bed mybedfile.bed --ref ref.fasta --mutect_args " --force_output --force_alleles " '
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --tumor_bam_folder   FOLDER                  Folder containing tumor BAM files to be called.'
    log.info '    --normal_bam_folder  FOLDER                  Folder containing matched normal BAM files.'
    log.info '    --bed                FILE                    Bed file containing intervals.'
    log.info '    --ref                FILE (with index)       Reference fasta file.'
    log.info '    --mutect_path        FILE                    mutect*.jar explicit path.'
    log.info '    --dbsnp              FILE                    dbSNP VCF file required by mutect.'
    log.info '    --cosmic             FILE                    Cosmic VCF file required by mutect.'
    log.info 'Optional arguments:'
    log.info '    --mutect_args        STRING                  Arguments you want to pass to mutect.'
    log.info '                                                 WARNING: form is " --force_alleles " with spaces between quotes.'  
    log.info '    --suffix_tumor       STRING                  Suffix identifying tumor bam (default: "_T").'
    log.info '    --suffix_normal      STRING                  Suffix identifying normal bam (default: "_N").'
    log.info '    --mem                INTEGER                 Java memory passed to mutect.'
    log.info '    --out_folder         FOLDER                  Output folder (default: mutect_results).'
    log.info ''
    log.info 'CAUTION: mutect*.jar executing file has to be in your $PATH.'
    log.info ''
    exit 1
}

fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )
fasta_ref_gzi = file( params.ref+'.gzi' )
fasta_ref_dict = file( params.ref.replace(".fasta",".dict") )
bed = file(params.bed)

params.suffix_tumor = "_T"
params.suffix_normal = "_N"
params.mem = 8
params.out_folder = "mutect_results"
params.mutect_args = ""

// FOR TUMOR 
// recovering of bam files
tumor_bams = Channel.fromPath( params.tumor_bam_folder+'/*'+params.suffix_tumor+'.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${params.tumor_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_tumor}.bam",""), path ] }

// recovering of bai files
tumor_bais = Channel.fromPath( params.tumor_bam_folder+'/*'+params.suffix_tumor+'.bam.bai' )
              .ifEmpty { error "Cannot find any bai file in: ${params.tumor_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_tumor}.bam.bai",""), path ] }

// building bam-bai pairs
tumor_bam_bai = tumor_bams
    	      .phase(tumor_bais)
    	      .map { tumor_bam, tumor_bai -> [ tumor_bam[0], tumor_bam[1], tumor_bai[1] ] }

// FOR NORMAL 
// recovering of bam files
normal_bams = Channel.fromPath( params.normal_bam_folder+'/*'+params.suffix_normal+'.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${params.normal_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_normal}.bam",""), path ] }

// recovering of bai files
normal_bais = Channel.fromPath( params.normal_bam_folder+'/*'+params.suffix_normal+'.bam.bai' )
              .ifEmpty { error "Cannot find any bai file in: ${params.normal_bam_folder}" }
              .map {  path -> [ path.name.replace("${params.suffix_normal}.bam.bai",""), path ] }

// building bam-bai pairs
normal_bam_bai = normal_bams
    	      .phase(normal_bais)
    	      .map { normal_bam, normal_bai -> [ normal_bam[0], normal_bam[1], normal_bai[1] ] }

// building 4-uplets corresponding to {tumor_bam, tumor_bai, normal_bam, normal_bai}
tn_bambai = tumor_bam_bai
	      .phase(normal_bam_bai)
	      .map {tumor_bb, normal_bb -> [ tumor_bb[1], tumor_bb[2], normal_bb[1], normal_bb[2] ] }	
// here each element X of tn_bambai channel is a 4-uplet. X[0] is the tumor bam, X[1] the tumor bai, X[2] the normal bam and X[3] the normal bai.

process mutect {
  
    tag { tumor_normal_tag }

    publishDir params.out_folder, mode: 'move'

    input:
    file tn from tn_bambai
    file bed
    file fasta_ref
    file fasta_ref_fai
    file fasta_ref_gzi
    file fasta_ref_dict

    output:
    file("${tumor_normal_tag}_calls.vcf") into mutect_output1
    file("${tumor_normal_tag}_calls_stats.txt") into mutect_output2

    shell:
    tumor_normal_tag = tn[0].baseName.replace(params.suffix_tumor,"")
    '''
    java -Xmx!{params.mem}g -jar !{params.mutect_path} --analysis_type MuTect --reference_sequence !{fasta_ref} --dbsnp !{params.dbsnp} --cosmic !{params.cosmic} --intervals !{bed} --input_file:tumor !{tumor_normal_tag}!{params.suffix_tumor}.bam --input_file:normal !{tumor_normal_tag}!{params.suffix_normal}.bam --out "!{tumor_normal_tag}_calls_stats.txt" --vcf "!{tumor_normal_tag}_calls.vcf" !{params.mutect_args}
    '''
}



