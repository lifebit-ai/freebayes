//============================================================================================================================
//        PARAMETERS
//============================================================================================================================


/*--------------------------------------------------
  Fasta related input files
  You can use the flag --hg19 for using the hg19 version of the Genome.
  You can use the flag --h38 for using the GRCh38.p10 version of the Genome.
  They can be passed manually, through the parameter:
  	params.fasta="/my/path/to/file";
  And if already at user's disposal ( if not, automatically generated ):
	params.fai="/my/path/to/file";
	params.fastagz="/my/path/to/file";
	params.gzfai="/my/path/to/file";
	params.gzi="/my/path/to/file";
  if no input is given, the hg19 version of the genome is used.
---------------------------------------------------*/
params.hg19="true";
params.h38="";
params.test="";
params.hg19chr20="";


params.fasta="nofasta";
params.fai="nofai";
params.fastagz="nofastagz";
params.gzfai="nogzfai";
params.gzi="nogzi";

if(!("nofasta").equals(params.fasta)){
  fasta=file(params.fasta)
  fai=file(params.fai);
  fastagz=file(params.fastagz);
  gzfai=file(params.gzfai);
  gzi=file(params.gzi);
}
else if(params.h38 ){
  fasta=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa");
  fai=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.fai");
  fastagz=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.gz");
  gzfai=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.gz.fai");
  gzi=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.gz.gzi");
}
else if(params.test ){
  fasta=file("s3://deepvariant-test/input/ucsc.hg19.chr20.unittest.fasta");
  fai=file("s3://deepvariant-test/input//ucsc.hg19.chr20.unittest.fasta.fai");
  fastagz=file("s3://deepvariant-test/input/ucsc.hg19.chr20.unittest.fasta.gz");
  gzfai=file("s3://deepvariant-test/input/ucsc.hg19.chr20.unittest.fasta.gz.fai");
  gzi=file("s3://deepvariant-test/input/ucsc.hg19.chr20.unittest.fasta.gz.gzi");
}
else if(params.hg19 ){
  fasta=file("s3://deepvariant-data/genomes/hg19/hg19.fa");
  fai=file("s3://deepvariant-data/genomes/hg19/hg19.fa.fai");
  fastagz=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz");
  gzfai=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz.fai");
  gzi=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz.gzi");
}
else if(params.hg19chr20 ){
  fasta=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa");
  fai=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa.fai");
  fastagz=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa.gz");
  gzfai=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa.gz.fai");
  gzi=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa.gz.gzi");
}
else{
  System.out.println(" --fasta \"/path/to/your/genome\"  params is required and was not found! ");
  System.out.println(" or you can use standard genome versions by typing --hg19 or --h38 ");
  System.exit(0);
}

int cores = Runtime.getRuntime().availableProcessors();
params.j=cores



/*--------------------------------------------------
  Params for the Read Group Line to be added just in
  case its needed.
  If not given, default values are used.
---------------------------------------------------*/
params.rgid=4;
params.rglb="lib1";
params.rgpl="illumina";
params.rgpu="unit1";
params.rgsm=20;

/*--------------------------------------------------
  Bam input files
  The input must be a path to a folder containing multiple bam files
---------------------------------------------------*/
params.bam_folder="s3://deepvariant-test/input/";
params.bam_file_prefix="*"


Channel.fromPath("${params.bam_folder}/${params.bam_file_prefix}*.bam").map{ file -> tuple(file.name, file) }.set{bamChannel}

/*--------------------------------------------------
  Variant Caller to be used
---------------------------------------------------*/
params.method="gatk-hc"

/*--------------------------------------------------
  Output directory
---------------------------------------------------*/
params.resultdir = "./Results";

//============================================================================================================================
//        PROCESSES
//============================================================================================================================

/******
*
*   PREPROCESSING:
*
*   1A Preprocessing Genome
*     - Index and compress the genome ( Create fa.fai, fa.gz, fa.gz.fai, fa.gz.gzi and the dict file )
*
*   1B Preprocessing Bams
*     - Add RG line in case it is missing
*     - Reorder Bam
*     - Index Bam
*
*
*****/

process preprocess_genome{

  container 'lifebitai/samtools'


  input:
  file fasta from fasta
  file fai from fai
  file fastagz from fastagz
  file gzfai from gzfai
  file gzi from gzi
  output:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"),file("${fasta.baseName}.dict") into fastaChannel
  script:
  """
  [[ ${params.fai} == "nofai" ]] &&  samtools faidx $fasta || echo " fai file of user is used, not created"
  [[ ${params.fastagz} == "nofastagz" ]]  && bgzip -c ${fasta} > ${fasta}.gz || echo "fasta.gz file of user is used, not created "
  [[ ${params.gzfai} == "nogzi" ]] && bgzip -c -i ${fasta} > ${fasta}.gz || echo "gzi file of user is used, not created"
  [[ ${params.gzi} == "nogzfai" ]] && samtools faidx "${fasta}.gz" || echo "gz.fai file of user is used, not created"
  java -jar /picard.jar CreateSequenceDictionary R= $fasta O= ${fasta.baseName}.dict
  """
}


process preprocess_bam {

  tag "${bam[0]}"
  container 'lifebitai/samtools'


  input:
  set val(prefix), file(bam) from bamChannel
  set file(genome),file(genomefai),file(genomegz),file(genomegzfai),file(genomegzgzi),file(genomedict) from fastaChannel

  output:
  set file("ready/${bam[0]}"), file("ready/${bam[0]}.bai") into completeChannel

  script:
  """
  ## Add RG line in case it is missing
    mkdir ready
    [[ `samtools view -H ${bam[0]} | grep '@RG' | wc -l`   > 0 ]] && { mv $bam ready; }|| { java -jar  /picard.jar  AddOrReplaceReadGroups \
    I=${bam[0]} O=ready/${bam[0]} RGID=${params.rgid} RGLB=${params.rglb} RGPL=${params.rgpl} RGPU=${params.rgpu} RGSM=${params.rgsm}; }
  ## Index Bam file
    cd ready; samtools index ${bam[0]};
  """

}


//Preparing channel ( pairing fasta with bams)
fastaChannel.map{file -> tuple (1,file[0],file[1],file[2],file[3],file[4],file[5])}
              .set{all_fa};

completeChannel.map { file -> tuple(1,file[0],file[1]) }
                 .set{all_bam};

all_fa.cross(all_bam)
        .set{all_fa_bam};

//all_fa_bam.subscribe{ println it };


/******
*
*   VARIANT CALLING
*
*****/

process run_variant_caller {

    tag "${bam[1]}"
    container "lifebitai/freebayes"
    publishDir "$baseDir/${params.resultdir}"
    cpus params.j

    input:
    set file(fasta), file(bam) from all_fa_bam

    output:
    file('calling_output.vcf') into methods_result

    script:
    """
    pwd=\$PWD
    cd ~/freebayes/scripts/
    ./freebayes-parallel <(./fasta_generate_regions.py \$pwd/${fasta[2]} 100000) ${task.cpus} -f \$pwd/${fasta[1]} \$pwd/${bam[1]} > \$pwd/calling_output.vcf
    """
}



workflow.onComplete {
    println ( workflow.success ? "Done! \nYou can find your results in $baseDir/${params.resultdir}" : "Oops .. something went wrong" )
}
