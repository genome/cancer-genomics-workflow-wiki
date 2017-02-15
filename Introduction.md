# Initialization

## Shell Variables

* NOTE: The directory where you download the following config files will become the base directory for this tutorial.
* You will need to update the $SOMATIC_HOME environment variable in Somatic_base_config.sh to match your local directory.

```bash

curl -L -k -o Somatic_config.sh https://raw.githubusercontent.com/genome/arvados_trial_wiki/master/Somatic_config.sh
curl -L -k -o Somatic_base_config.sh https://raw.githubusercontent.com/genome/arvados_trial_wiki/master/Somatic_base_config.sh

```

* NOTE: This must be run each time a new terminal window is opened.  These shell variables will NOT persist in your environment.

```bash

source Somatic_config.sh

```
## Working Directories

```bash

mkdir -p $SOMATIC_HOME/software
mkdir -p $SOMATIC_HOME/logs
mkdir -p $SOMATIC_HOME/refseq
mkdir -p $SOMATIC_HOME/fastq
mkdir -p $SOMATIC_HOME/alignments
mkdir -p $SOMATIC_HOME/varscan
mkdir -p $SOMATIC_HOME/strelka/exome
mkdir -p $SOMATIC_HOME/strelka/wgs
mkdir -p $SOMATIC_HOME/mutect/exome
mkdir -p $SOMATIC_HOME/mutect/wgs
mkdir -p $SOMATIC_HOME/pindel
mkdir -p $VEP_CACHE

```

# Installs and Downloads

## Software

### [gsutil](https://cloud.google.com/storage/docs/gsutil_install)

This tool is required to download files from the Broad Google Cloud Platform
For install instructions for your specific architecture, please see the above link to gsutil.


NOTE: THIS DID NOT WORK FOR ME ON LINUX WORKSTATION, I INSTALLED ON MY MAC

```bash

sudo easy_install -U pip
sudo pip install gsutil

```

### [HTSlib](https://github.com/samtools/htslib)

```bash

cd $SOMATIC_HOME/software

curl -L -k -o htslib-1.3.2.tar.bz2 https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2
tar --bzip2 -xvf htslib-1.3.2.tar.bz2
cd htslib-1.3.2
./configure  --enable-plugins --enable-libcurl --prefix=$SOMATIC_HOME/software
make
make install

$SOMATIC_HOME/software/bin/tabix 

```

### [Samtools](https://github.com/samtools/samtools)

```bash

cd $SOMATIC_HOME/software
curl -L -k -o samtools-1.3.1.tar.bz2  https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar --bzip2 -xvf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
./configure --with-htslib=$SOMATIC_HOME/software/htslib-1.3.2 --prefix=$SOMATIC_HOME/software
make
make install

$SOMATIC_HOME/software/bin/samtools

```

### [BWA](https://sourceforge.net/projects/bio-bwa/)

```bash

cd $SOMATIC_HOME/software
curl -L -k -o bwa-0.7.15.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2/download
tar --bzip2 -xvf bwa-0.7.15.tar.bz2
cd  bwa-0.7.15
make

ln -s $SOMATIC_HOME/software/bwa-0.7.15/bwa $SOMATIC_HOME/software/bin/bwa

bwa

```

### [Samblaster](https://github.com/GregoryFaust/samblaster)

```bash

cd $SOMATIC_HOME/software
curl -L -k -o samblaster-v.0.1.24.tar.gz https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.24/samblaster-v.0.1.24.tar.gz
tar -xzvf samblaster-v.0.1.24.tar.gz
cd samblaster-v.0.1.24
make

ln -s $SOMATIC_HOME/software/samblaster-v.0.1.24/samblaster $SOMATIC_HOME/software/bin/samblaster

samblaster -h

```

### [Sambamba](http://lomereiter.github.io/sambamba/)

```bash

cd $SOMATIC_HOME/software
curl -L -k -o sambamba_v0.6.4_linux.tar.bz2 https://github.com/lomereiter/sambamba/releases/download/v0.6.4/sambamba_v0.6.4_linux.tar.bz2
tar --bzip2 -xvf sambamba_v0.6.4_linux.tar.bz2

ln -s $SOMATIC_HOME/software/sambamba_v0.6.4 $SOMATIC_HOME/software/bin/sambamba

sambamba -h

```

### [Picard](https://broadinstitute.github.io/picard/)

```bash

cd $SOMATIC_HOME/software
curl -L -k -o picard-tools-2.4.1.zip https://github.com/broadinstitute/picard/releases/download/2.4.1/picard-tools-2.4.1.zip
unzip picard-tools-2.4.1.zip

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar -h 

```

### [GATK](https://software.broadinstitute.org/gatk/)

1. Manually download GATK after accepting the license : https://www.broadinstitute.org/gatk/download/auth?package=GATK
2. Copy the download : `scp Downloads/GenomeAnalysisTK-3.6.tar.bz2 linus2112:$SOMATIC_HOME/software/.`
3. Unzip the archive : `tar --bzip2 -xvf GenomeAnalysisTK-3.6.tar.bz2`
4. Test with JAVA8 : `$JAVA_EIGHT -jar GenomeAnalysisTK.jar -h`

### [verifyBamId](http://genome.sph.umich.edu/wiki/VerifyBamID)

```bash

cd $SOMATIC_HOME/software
curl -L -k -o verifyBamIDLibStatGen.1.1.3.tgz https://github.com/statgen/verifyBamID/releases/download/v1.1.3/verifyBamIDLibStatGen.1.1.3.tgz
tar -xzvf verifyBamIDLibStatGen.1.1.3.tgz
cd verifyBamID_1.1.3
make

ln -s $SOMATIC_HOME/software/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID $SOMATIC_HOME/software/bin/verifyBamID

verifyBamID

```

### [VarScan](http://dkoboldt.github.io/varscan/)

```bash

cd $SOMATIC_HOME/software
curl -L -k -o VarScan.v2.4.2.jar https://github.com/dkoboldt/varscan/releases/download/2.4.2/VarScan.v2.4.2.jar
$JAVA_EIGHT -jar $SOMATIC_HOME/software/VarScan.v2.4.2.jar

```

### [Strelka](https://sites.google.com/site/strelkasomaticvariantcaller/)

```bash

cd $SOMATIC_HOME/software

curl -L -k -o strelka-2.7.1.centos5_x86_64.tar.bz2 https://github.com/Illumina/strelka/releases/download/v2.7.1/strelka-2.7.1.centos5_x86_64.tar.bz2
tar --bzip2 -xvf strelka-2.7.1.centos5_x86_64.tar.bz2
bash $SOMATIC_HOME/software/strelka-2.7.1.centos5_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash

```

### [BCFtools](https://samtools.github.io/bcftools/)

```bash

cd $SOMATIC_HOME/software
curl -L -k -o bcftools-1.3.1.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
tar --bzip2 -xvf bcftools-1.3.1.tar.bz2
cd bcftools-1.3.1
make prefix=$SOMATIC_HOME/software install
cd ../
bcftools -h

```

### [Pindel](https://github.com/genome/pindel)

```bash

cd $SOMATIC_HOME/software
curl -L -k -o pindel-0.2.5b8.tar.gz https://github.com/genome/pindel/archive/v0.2.5b8.tar.gz
tar -xzvf pindel-0.2.5b8.tar.gz
cd pindel-0.2.5b8/
./INSTALL /gscuser/jwalker/git/HCC1395/arvados/software/htslib-1.3.2
ln -s $SOMATIC_HOME/software/pindel-0.2.5b8/pindel $SOMATIC_HOME/software/bin/pindel
ln -s $SOMATIC_HOME/software/pindel-0.2.5b8/pindel2vcf $SOMATIC_HOME/software/bin/pindel2vcf
cd $SOMATIC_HOME
pindel -h

```

### [VEP](https://github.com/Ensembl/ensembl-vep)

```bash

cd $SOMATIC_HOME/software/
git clone https://github.com/Ensembl/ensembl-vep.git --branch release/87 --single-branch
cd $SOMATIC_HOME/software/ensembl-vep
unset PERL5LIB
/usr/bin/perl INSTALL.pl --NO_HTSLIB --CACHEDIR $VEP_CACHE
#During the installation make sure to accept (y) when asked whether you'd like to install cache files, fastas, and plugins.
#Install all homo sapiens build 38 cache files (options 42 44 46), the homo sapiens fasta (option 28), and the Downstream plugin (option 24).
wget -O $VEP_CACHE/Plugins/Wildtype.pm https://raw.githubusercontent.com/griffithlab/pVAC-Seq/master/pvacseq/VEP_plugins/Wildtype.pm --no-check-certificate


```

### [bam-readcount](https://github.com/genome/bam-readcount)

* NOTE: This requires cmake 2.8.3 or higher

```bash

cd $SOMATIC_HOME/software/
mkdir git
cd git
git clone --recursive git://github.com/genome/bam-readcount.git

cd $SOMATIC_HOME/software/
mkdir bam-readcount
cd bam-readcount
export SAMTOOLS_ROOT=$SOMATIC_HOME/software/samtools-1.3.1
cmake $SOMATIC_HOME/software/git/bam-readcount
make
./bin/bam-readcount

```

### [fpfilter](https://github.com/genome/fpfilter-tool)

```bash

cd $SOMATIC_HOME/software/
curl -L -k -o fpfilter.pl https://raw.githubusercontent.com/genome/fpfilter-tool/master/fpfilter.pl
/usr/bin/perl $SOMATIC_HOME/software/fpfilter.pl

```

## Reference Sequence

### [1000 Genomes](http://www.internationalgenome.org/)

```bash

mkdir -p $SOMATIC_REFSEQ_DIR

cd $SOMATIC_$REFSEQ_DIR

curl -L -k -o $SOMATIC_REFSEQ_FASTA  $GENOME_URI

```

### [GATK Resources](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/?pli=1)

TODO: Make WGS intervals using Picard tools

```bash 

cd $SOMATIC_REFSEQ_DIR

gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/$REFSEQ_DBSNP $SOMATIC_REFSEQ_DIR
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/$REFSEQ_KNOWN_INDELS $SOMATIC_REFSEQ_DIR
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/$REFSEQ_MILLS_INDELS $SOMATIC_REFSEQ_DIR

gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list $SOMATIC_REFSEQ_DIR
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ $SOMATIC_REFSEQ_DIR

```

#### [COSMIC](http://cancer.sanger.ac.uk/cosmic)

```bash

cd $SOMATIC_REFSEQ_DIR

sftp "your_email_address"@sftp-cancer.sanger.ac.uk
get /files/grch38/cosmic/v79/VCF/CosmicCodingMuts.vcf.gz
get /files/grch38/cosmic/v79/VCF/CosmicNonCodingVariants.vcf.gz

zgrep "^#" CosmicCodingMuts.vcf.gz > VCF_Header
zgrep -v "^#" CosmicCodingMuts.vcf.gz | awk '{print "chr"$0}' | sed 's/^chrMT/chrM/' > CosmicCodingMuts.clean
zgrep -v "^#" CosmicNonCodingVariants.vcf.gz | awk '{print "chr"$0}' | sed 's/^chrMT/chrM/' > CosmicNonCodingVariants.clean

cat CosmicCodingMuts.clean CosmicNonCodingVariants.clean | sort -gk 2,2 > Cosmic_v79
cat VCF_Header Cosmic_v79 > Cosmic_v79.vcf

rm VCF_Header CosmicCodingMuts.clean CosmicNonCodingVariants.clean Cosmic_v79

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar SortVcf I=$SOMATIC_REFSEQ_DIR/Cosmic_v79.vcf O=$SOMATIC_REFSEQ_DIR/Cosmic_v79.dictsorted.vcf SEQUENCE_DICTIONARY=$SOMATIC_REFSEQ_DICT

```

### Exome Targets/Probes

#### [IDT xGen Exome Research Panel v1.0](http://www.idtdna.com/pages/products/nextgen/target-capture/xgen-lockdown-panels/xgen-exome-panel)

```bash

cd $SOMATIC_HOME/refseq/

curl -L -k -o xgen-exome-research-panel-probes.bed http://www.idtdna.com/pages/docs/default-source/xgen-libraries/xGen-Lockdown-Panels/xgen-exome-research-panel-probes.bed?sfvrsn=4

curl -L -k -o xgen-exome-research-panel-targets.bed http://www.idtdna.com/pages/docs/default-source/xgen-libraries/xGen-Lockdown-Panels/xgen-exome-research-panel-targets.bed?sfvrsn=6

```

#### [UCSC Genome Browser](https://genome.ucsc.edu/)
Download the UCSC hg19ToHg38 chain file

```base

cd $SOMATIC_HOME/refseq/

curl -L -k -o hg19ToHg38.over.chain.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz

```

# Index

## Reference Sequence

### BWA Alinger Index

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/bwa_index.out -e $SOMATIC_HOME/logs/bwa_index.err -M ${INDEX_RAM_GB}000000 -R "select[mem>=${INDEX_RAM_GB}000] rusage[mem=${INDEX_RAM_GB}000]" $SOMATIC_HOME/software/bwa-0.7.15/bwa index $SOMATIC_REFSEQ_FASTA

```

### Picard Sequence Dictionary

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/picard_dict.out -e $SOMATIC_HOME/logs/picard_dict.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=$SOMATIC_REFSEQ_FASTA O=$SOMATIC_REFSEQ_DICT GENOME_ASSEMBLY=$GENOME_ASSEMBLY URI=$GENOME_URI SPECIES=$GENOME_SPECIES

```

### Samtools FAIDX

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/samtools_faidx.out -e $SOMATIC_HOME/logs/samtools_faidx.err $SOMATIC_HOME/software/bin/samtools faidx $SOMATIC_REFSEQ_FASTA

```

## GATK Known Sites

### [Tabix](http://www.htslib.org/doc/tabix.html)

#### dbSNP VCF

```bash

bgzip $SOMATIC_REFSEQ_DIR/Homo_sapiens_assembly38.dbsnp138.vcf
tabix -p vcf $SOMATIC_REFSEQ_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.gz

```

#### Known Indels

```bash

tabix -p vcf $SOMATIC_REFSEQ_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz

```

#### Mills Indels

```bash

tabix -p vcf $SOMATIC_REFSEQ_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

```

#### COSMIC VCF

```bash

bgzip $SOMATIC_REFSEQ_DIR/Cosmic_v79.dictsorted.vcf
tabix -p vcf $SOMATIC_REFSEQ_DIR/Cosmic_v79.dictsorted.vcf.gz

```

## Target Interval Lists

### Exome

NOTE: hg19 sequence dictionary used from GMS @ MGI

### [Picard BedToIntervalList](https://broadinstitute.github.io/picard/command-line-overview.html#BedToIntervalList)

```bash

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar BedToIntervalList I=$SOMATIC_HOME/refseq/xgen-exome-research-panel-probes.bed SD=/gscmnt/gc4096/info/model_data/2871743894/build108563338/seqdict/seqdict.sam O=$SOMATIC_HOME/refseq/xgen-exome-research-panel-probes.interval_list

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar BedToIntervalList I=$SOMATIC_HOME/refseq/xgen-exome-research-panel-targets.bed SD=/gscmnt/gc4096/info/model_data/2871743894/build108563338/seqdict/seqdict.sam O=$SOMATIC_HOME/refseq/xgen-exome-research-panel-targets.interval_list

```

### [Picard LiftOverIntervalList](https://broadinstitute.github.io/picard/command-line-overview.html#LiftOverIntervalList)

```bash

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar LiftOverIntervalList CHAIN=$SOMATIC_HOME/refseq/hg19ToHg38.over.chain SD=$SOMATIC_REFSEQ_DICT O=$SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-probes.interval_list I=$SOMATIC_HOME/refseq/xgen-exome-research-panel-probes.interval_list

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar LiftOverIntervalList CHAIN=$SOMATIC_HOME/refseq/hg19ToHg38.over.chain SD=$SOMATIC_REFSEQ_DICT O=$SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-targets.interval_list I=$SOMATIC_HOME/refseq/xgen-exome-research-panel-targets.interval_list

```

### IntervalListToBed using [Perl](https://www.perl.org/)

```bash

cat $SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-probes.interval_list | /usr/bin/perl -M5.10.0 -ne 'if(substr($_,0,1) ne q{@}) { chomp; my @c = split "\t"; say(join("\t", $c[0], $c[1]-1, $c[2])); }' > $SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-probes.bed

cat $SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-targets.interval_list | /usr/bin/perl -M5.10.0 -ne 'if(substr($_,0,1) ne q{@}) { chomp; my @c = split "\t"; say(join("\t", $c[0], $c[1]-1, $c[2])); }' > $SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-targets.bed

```

### Scatter with [Picard IntervalListTools](https://broadinstitute.github.io/picard/command-line-overview.html#IntervalListTools)

```bash

mkdir -p $SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-targets

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar IntervalListTools INPUT=$SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-targets.interval_list SCATTER_COUNT=$REFERENCE_SEQUENCE_CHUNKS OUTPUT=$SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-targets

```

### WGS

#### Whole-Genome Interval List


```bash

awk '{print $1"\t1\t"$2"\t+\t"$1}' ${SOMATIC_REFSEQ_FASTA}.fai | cat $SOMATIC_REFSEQ_DICT - > $SOMATIC_REFSEQ_DIR/$GENOME_BASENAME.interval_list
)

```

#### Autosomal Chromosome Interval List

```bash

egrep 'chr[0-9]{1,2}\s' ${SOMATIC_REFSEQ_FASTA}.fai | awk '{print $1"\t1\t"$2"\t+\t"$1}' | cat $SOMATIC_REFSEQ_DICT - > $SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_autosomal.interval_list

```

#### Autosomal + Sex Chromosome Interval List

```bash

egrep 'chr[0-9,X,Y]{1,2}\s' ${SOMATIC_REFSEQ_FASTA}.fai | awk '{print $1"\t1\t"$2"\t+\t"$1}' | cat $SOMATIC_REFSEQ_DICT - > $SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_autosomal_plus_sex.interval_list

```

#### [Picard ScatterIntervalsByNs](https://broadinstitute.github.io/picard/command-line-overview.html#ScatterIntervalsByNs)

NOTE : It looks like Broad may use a higher value for Ns (>1) when making a list of calling regions. There are only 356 regions compared to 592 when using 1. 100 for MAX_TO_MERGE resulted in 324 regions and 99 resulted in 533??? Try kdiff3 tomorrow.

```base

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar ScatterIntervalsByNs R=$SOMATIC_REFSEQ_FASTA O=$SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_excludeNs.interval_list OT=ACGT N=100

```

#### [Picard IntervalListTools](https://broadinstitute.github.io/picard/command-line-overview.html#IntervalListTools) Intersect ExcludeNs with Autosomal + Sex Interval List

```bash

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar IntervalListTools I=$SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_excludeNs.interval_list I=$SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_autosomal_plus_sex.interval_list ACTION=INTERSECT O=$SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_calling_regions.interval_list

```

```bash

$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar IntervalListTools I=$SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_calling_regions.interval_list SCATTER_COUNT=$REFERENCE_SEQUENCE_CHUNKS OUTPUT=$SOMATIC_REFSEQ_DIR

```


