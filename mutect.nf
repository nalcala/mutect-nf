#! /usr/bin/env nextflow

params.help = null

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
    log.info '    OR'
    log.info '    --bam_folder         FILE                    Folder containing all bam files.'
    log.info '    --tn_file            FILE                    Tab delimited text file with two columns called tumor and normal'
    log.info '                                                 where each line contains the path of two matched BAM files.'    
    log.info '    --ref                FILE (with index)       Reference fasta file.'
    log.info '    --mutect_jar         FILE                    mutect*.jar explicit path.'
    log.info '    OR'
    log.info '    --mutect2_jar        FILE                    gatk*.jar explicit path to run mutect2 (integrated to GATK).'
    log.info 'Optional arguments:'
    log.info '    --bed                FILE                    Bed file containing intervals.'
    log.info '    --region             REGION                  A region defining the calling, in the format CHR:START-END.'
    log.info '    NOTE: if neither --bed or --region, will perform the calling on whole genome, based on the faidx file.'
    log.info '    --nsplit             INTEGER                 Split the region for calling in nsplit pieces and run in parallel.'
    log.info '    --dbsnp              FILE                    dbSNP VCF file required by mutect (CAUTION: not a symbolic link).'
    log.info '    --cosmic             FILE                    Cosmic VCF file required by mutect (CAUTION: not a symbolic link).'
    log.info '    --mutect_args        STRING                  Arguments you want to pass to mutect.'
    log.info '                                                 WARNING: form is " --force_alleles " with spaces between quotes.'
    log.info '    --suffix_tumor       STRING                  Suffix identifying tumor bam (default: "_T").'
    log.info '    --suffix_normal      STRING                  Suffix identifying normal bam (default: "_N").'
    log.info '    --mem                INTEGER                 Java memory passed to mutect.'
    log.info '    --out_folder         FOLDER                  Output folder (default: mutect_results).'
    log.info '    --java               PATH                    Name of the JAVA command  (default: java).'
    log.info ''
    log.info ''
    exit 0
}

fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )
fasta_ref_gzi = file( params.ref+'.gzi' )
fasta_ref_dict = file( params.ref.replace(".fasta",".dict").replace(".fa",".dict") )

params.suffix_tumor = "_T"
params.suffix_normal = "_N"
params.mem = 8
params.out_folder = "mutect_results"
params.mutect_args = ""
params.nsplit = 1
params.region = null
params.bed = null
params.java = "java"
params.dbsnp = ""
params.cosmic = ""
if (params.dbsnp == "") { dbsnp_option = "" } else { dbsnp_option = "--dbsnp" }
if (params.cosmic == "") { cosmic_option = "" } else { cosmic_option = "--cosmic" }
params.mutect_jar = null
params.mutect2_jar = null
mutect_version = params.mutect_jar ? 1 : 2
params.tn_file = null
params.bam_folder = "./"
params.tumor_bam_folder = null
params.normal_bam_folder = null

if (params.tn_file) {
    // FOR INPUT AS A TAB DELIMITED FILE
    tn_bambai = Channel.fromPath(params.tn_file).splitCsv(header: true, sep: '\t', strip: true).map{row -> [ file(params.bam_folder + "/" + row.tumor), file(params.bam_folder + "/" + row.tumor+'.bai') ,file(params.bam_folder + "/" + row.normal), file(params.bam_folder + "/" + row.normal+'.bai') ]}
} else {
    // FOR INPUT AS TWO FOLDER
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
}

/* manage input positions to call (bed or region or whole-genome) */
  if (params.region) {
      input_region = 'region'
  } else if (params.bed) {
      input_region = 'bed'
      bed = file(params.bed)
  } else {
      input_region = 'whole_genome'
}

/* process to create a bed file from region or from faidx if whole-genome, otherwise return the input bed file */
process bed {
      output:
      file "temp.bed" into outbed

      shell:
      if (input_region == 'region')
      '''
      echo !{params.region} | sed -e 's/[:|-]/	/g' > temp.bed
      '''

      else if (input_region == 'bed')
      '''
      ln -s !{bed} temp.bed
      '''

      else if (input_region == 'whole_genome')
      '''
      cat !{fasta_ref_fai} | awk '{print $1"	"1"	"$2 }' > temp.bed
      '''
  }


/* split bed file into nsplit regions */
process split_bed {

      input:
      file bed from outbed

      output:
      file '*_regions.bed' into split_bed, count_split_bed mode flatten

      shell:
      '''
      grep -v '^track' !{bed} | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print $1" "$2" "$3}' | cut_into_small_beds.r !{params.nsplit}
      '''
}

//println count_split_bed.count().val

process mutect {

    memory params.mem+'GB' 

    tag { printed_tag }

    input:
    val mutect_version
    val dbsnp_option
    val cosmic_option
    file bed_tn from tn_bambai.spread(split_bed)
    file fasta_ref
    file fasta_ref_fai
    file fasta_ref_gzi
    file fasta_ref_dict

    output:
    set val(tumor_normal_tag), file("${tumor_normal_tag}_${bed_tag}_calls.vcf") into mutect_output1
    set val(tumor_normal_tag), file("${tumor_normal_tag}_${bed_tag}_calls_stats.txt") into mutect_output2

    shell:
    tumor_normal_tag = bed_tn[0].baseName.replace(params.suffix_tumor,"")
    tumor_normal_tag = bed_tn[0].baseName.replace(params.suffix_tumor,"") + "_vs_" + bed_tn[2].baseName.replace(params.suffix_normal,"")
    bed_tag = bed_tn[4].baseName //bed_tn = bamN,baiN,bamT,baiT,bed
    printed_tag = tumor_normal_tag + "_" + bed_tag
    '''
    if [ "!{mutect_version}" == "2" ]
        then
            !{params.java} -Xmx!{params.mem}g -jar !{params.mutect2_jar} -T MuTect2 -R !{fasta_ref} !{dbsnp_option} !{params.dbsnp} !{cosmic_option} !{params.cosmic} -I:normal !{bed_tn[2]} -I:tumor !{bed_tn[0]} -o "!{tumor_normal_tag}_!{bed_tag}_calls.vcf" -L !{bed_tn[4]} !{params.mutect_args}
	    touch !{tumor_normal_tag}_!{bed_tag}_calls_stats.txt 
        else
            !{params.java} -Xmx!{params.mem}g -jar !{params.mutect_jar} --analysis_type MuTect --reference_sequence !{fasta_ref} !{dbsnp_option} !{params.dbsnp} !{cosmic_option} !{params.cosmic} --intervals !{bed_tn[4]} --input_file:tumor !{bed_tn[0]} --input_file:normal !{bed_tn[2]} --out "!{tumor_normal_tag}_!{bed_tag}_calls_stats.txt" --vcf "!{tumor_normal_tag}_!{bed_tag}_calls.vcf" !{params.mutect_args}
    fi
    '''
}

beds_length = count_split_bed.count().val

process mergeMuTectOutputs {

    tag { tumor_normal_tag }

    publishDir params.out_folder, mode: 'move'

    input:
    set val(tumor_normal_tag), file(vcf_files) from mutect_output1.groupTuple(size: beds_length)
    set val(tumor_normal_tag2), file(txt_files) from mutect_output2.groupTuple(size: beds_length)

    output:
    file("${tumor_normal_tag}_calls.vcf") into res1
    file("${tumor_normal_tag}_calls_stats.txt") optional true into res2

    shell:
    '''
    # MERGE VCF FILES
    sed '/^#CHROM/q' `ls -1 *.vcf | head -1` > header.txt
    # Check if sort command allows sorting in natural order (chr1 chr2 chr10 instead of chr1 chr10 chr2)
    if [ `sort --help | grep -c 'version-sort' ` == 0 ]
    then
        sort_ops="-k1,1d"
    else
        sort_ops="-k1,1V"
    fi
    # Add all VCF contents and sort
    grep --no-filename -v '^#' *.vcf | LC_ALL=C sort -t '	' $sort_ops -k2,2n >> header.txt
    mv header.txt !{tumor_normal_tag}_calls.vcf

    # MERGE TXT FILES if mutect1
    if [ "!{mutect_version}" != "2" ]
        then
          head -n2 `ls -1 *.txt | head -1` > header.txt
          sed -i '1,2d' *calls_stats.txt
          cat *calls_stats.txt >> header.txt
          mv header.txt !{tumor_normal_tag}_calls_stats.txt
    fi  
    '''
}
