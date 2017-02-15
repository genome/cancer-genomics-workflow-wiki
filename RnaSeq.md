# Initialize

## Environment and Shell Variables

```bash

source /path/to/RnaSeq_config.sh

```

## Make Directories

```bash

mkdir -p $RNA_HOME/software/src
mkdir -p $RNA_HOME/software/bin
mkdir -p $RNA_HOME/annotation
mkdir -p $RNA_HOME/refseq

mkdir -p $RNA_HOME/logs
mkdir -p $RNA_HOME/alignments
mkdir -p $RNA_HOME/transcripts
mkdir -p $RNA_HOME/abundance

```

## Software Install

Time ~30 minutes

### HISAT

```bash

cd $RNA_HOME/software
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
unzip hisat2-2.0.4-Linux_x86_64.zip
ln -s $RNA_HOME/software/hisat2-2.0.4/hisat2 $RNA_HOME/software/bin/hisat2
hisat2

```

### Sambamba

```bash

cd $RNA_HOME/software
curl -L -k -o sambamba_v0.6.4_linux.tar.bz2 https://github.com/lomereiter/sambamba/releases/download/v0.6.4/sambamba_v0.6.4_linux.tar.bz2
tar --bzip2 -xvf sambamba_v0.6.4_linux.tar.bz2
ln -s $RNA_HOME/software/sambamba_v0.6.4 $RNA_HOME/software/sambamba
./sambamba

```

### StringTie

```bash

cd $RNA_HOME/software
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz
tar -xzvf stringtie-1.3.0.Linux_x86_64.tar.gz
ln -s $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie $RNA_HOME/software/bin/stringtie
stringtie

```

### gffcompare

```bash

cd $RNA_HOME/software
wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.9.8.Linux_x86_64.tar.gz
tar -xzvf gffcompare-0.9.8.Linux_x86_64.tar.gz
ln -s $RNA_HOME/software/gffcompare-0.9.8.Linux_x86_64/gffcompare $RNA_HOME/software/bin/gffcompare 
gffcompare

```

### R

```bash

cd $RNA_HOME/software
export R_LIBS=
wget https://cran.r-project.org/src/base/R-3/R-3.2.5.tar.gz
tar -zxvf R-3.2.5.tar.gz
cd R-3.2.5
./configure --prefix=$RNA_HOME/software --with-x=no
make
make install
Rscript --version

```

#### Ballgown

First start an interactive R session

```bash

$RNA_HOME/software/bin/R

```

The following commands are executed from within the R session. 

```bash

install.packages("devtools",repos="http://cran.us.r-project.org")
source("http://bioconductor.org/biocLite.R")
biocLite(c("alyssafrazee/RSkittleBrewer","dplyr","genefilter","ballgown"))
quit()

```

# Setup

## Convert Unaligned BAM to FASTQ

Time ~75 minutes

```bash

bsub -o $RNA_HOME/logs/picard_sam-to-fastq_${TUMOR_DATA_1_BASE}.out -e $RNA_HOME/logs/picard_sam-to-fastq_${TUMOR_DATA_1_BASE}.err gmt picard standard-sam-to-fastq --use-version=1.123 --input=$TUMOR_DATA_1_BAM --fastq ${RNA_HOME}/fastq/${TUMOR_DATA_1_BASE}_1.fastq --second-end-fastq ${RNA_HOME}/fastq/${TUMOR_DATA_1_BASE}_2.fastq

bsub -o $RNA_HOME/logs/picard_sam-to-fastq_${TUMOR_DATA_2_BASE}.out -e $RNA_HOME/logs/picard_sam-to-fastq_${TUMOR_DATA_2_BASE}.err gmt picard standard-sam-to-fastq --use-version=1.123 --input=$TUMOR_DATA_2_BAM --fastq ${RNA_HOME}/fastq/${TUMOR_DATA_2_BASE}_1.fastq --second-end-fastq ${RNA_HOME}/fastq/${TUMOR_DATA_2_BASE}_2.fastq

bsub -o $RNA_HOME/logs/picard_sam-to-fastq_${NORMAL_DATA_1_BASE}.out -e $RNA_HOME/logs/picard_sam-to-fastq_${NORMAL_DATA_1_BASE}.err gmt picard standard-sam-to-fastq --use-version=1.123 --input=$NORMAL_DATA_1_BAM --fastq ${RNA_HOME}/fastq/${NORMAL_DATA_1_BASE}_1.fastq --second-end-fastq ${RNA_HOME}/fastq/${NORMAL_DATA_1_BASE}_2.fastq

bsub -o $RNA_HOME/logs/picard_sam-to-fastq_${NORMAL_DATA_2_BASE}.out -e $RNA_HOME/logs/picard_sam-to-fastq_${NORMAL_DATA_2_BASE}.err gmt picard standard-sam-to-fastq --use-version=1.123 --input=$NORMAL_DATA_2_BAM --fastq ${RNA_HOME}/fastq/${NORMAL_DATA_2_BASE}_1.fastq --second-end-fastq ${RNA_HOME}/fastq/${NORMAL_DATA_2_BASE}_2.fastq

```

## Annotation

Time ~10 hours (~1 hour on production apipe-builder/lims-pd-long)

### Transcriptome

NOTE : The following python commands will not work with /usr/bin/python on MGI hosts

```bash

$RNA_HOME/software/hisat2-2.0.4/hisat2_extract_splice_sites.py $REFERENCE_GTF > $RNA_HOME/annotation/${REFERENCE_NAME}_ss.tsv
$RNA_HOME/software/hisat2-2.0.4/hisat2_extract_exons.py $REFERENCE_GTF > $RNA_HOME/annotation/${REFERENCE_NAME}_exons.tsv

```

### Index

```bash

bsub -q bigmem -o $RNA_HOME/logs/hisat2_build.out -e $RNA_HOME/logs/hisat2_build.err -M ${INDEX_RAM_GB}000000 -R "select[mem>=${INDEX_RAM_GB}000] rusage[mem=${INDEX_RAM_GB}000] span[hosts=1]" -n $THREADS $RNA_HOME/software/hisat2-2.0.4/hisat2-build -p $THREADS --ss $RNA_HOME/annotation/${REFERENCE_NAME}_ss.tsv --exon $RNA_HOME/annotation/${REFERENCE_NAME}_exons.tsv $REFERENCE_FASTA $RNA_HOME/refseq/${REFERENCE_NAME}_tran

```

# RNA-seq Workflow

## Alignment

```bash

# Make TEMP directories
TUMOR_DATA_1_TEMP=`mktemp -d $RNA_HOME/alignments/${TUMOR_DATA_1_BASE}.XXXXXXXXXXXX`
TUMOR_DATA_2_TEMP=`mktemp -d $RNA_HOME/alignments/${TUMOR_DATA_2_BASE}.XXXXXXXXXXXX`
NORMAL_DATA_1_TEMP=`mktemp -d $RNA_HOME/alignments/${NORMAL_DATA_1_BASE}.XXXXXXXXXXXX`
NORMAL_DATA_2_TEMP=`mktemp -d $RNA_HOME/alignments/${NORMAL_DATA_2_BASE}.XXXXXXXXXXXX`

bsub -o $RNA_HOME/logs/hisat2_align_${TUMOR_DATA_1_BASE}.out -e $RNA_HOME/logs/hisat2_align_${TUMOR_DATA_1_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n $THREADS "$RNA_HOME/software/hisat2-2.0.4/hisat2 -p $THREADS --dta -x $RNA_HOME/refseq/${REFERENCE_NAME}_tran --rg-id $TUMOR_DATA_1_ID --rg PL:ILLUMINA --rg PU:${TUMOR_DATA_1_FC}-${TUMOR_DATA_1_BC}.${TUMOR_DATA_1_LN} --rg LB:${TUMOR_DATA_1_LB} --rg SM:${TUMOR_DATA_SM} --rna-strandness RF -1 $RNA_HOME/fastq/${TUMOR_DATA_1_BASE}_1.fastq -2 $RNA_HOME/fastq/${TUMOR_DATA_1_BASE}_2.fastq | $RNA_HOME/software/sambamba_v0.6.4 view -S -f bam -l 0 /dev/stdin | $RNA_HOME/software/sambamba_v0.6.4 sort -t $THREADS -m ${SAMBAMBA_RAM_GB}G --tmpdir $TUMOR_DATA_1_TEMP -o $RNA_HOME/alignments/${TUMOR_DATA_1_BASE}.bam /dev/stdin"

bsub -o $RNA_HOME/logs/hisat2_align_${TUMOR_DATA_2_BASE}.out -e $RNA_HOME/logs/hisat2_align_${TUMOR_DATA_2_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n $THREADS "$RNA_HOME/software/hisat2-2.0.4/hisat2 -p $THREADS --dta -x $RNA_HOME/refseq/${REFERENCE_NAME}_tran --rg-id $TUMOR_DATA_2_ID --rg PL:ILLUMINA --rg PU:${TUMOR_DATA_2_FC}-${TUMOR_DATA_2_BC}.${TUMOR_DATA_2_LN} --rg LB:${TUMOR_DATA_2_LB} --rg SM:${TUMOR_DATA_SM} --rna-strandness RF -1 $RNA_HOME/fastq/${TUMOR_DATA_2_BASE}_1.fastq -2 $RNA_HOME/fastq/${TUMOR_DATA_2_BASE}_2.fastq | $RNA_HOME/software/sambamba_v0.6.4 view -S -f bam -l 0 /dev/stdin | $RNA_HOME/software/sambamba_v0.6.4 sort -t $THREADS -m ${SAMBAMBA_RAM_GB}G --tmpdir $TUMOR_DATA_2_TEMP -o $RNA_HOME/alignments/${TUMOR_DATA_2_BASE}.bam /dev/stdin"

bsub -o $RNA_HOME/logs/hisat2_align_${NORMAL_DATA_1_BASE}.out -e $RNA_HOME/logs/hisat2_align_${NORMAL_DATA_1_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n $THREADS "$RNA_HOME/software/hisat2-2.0.4/hisat2 -p $THREADS --dta -x $RNA_HOME/refseq/${REFERENCE_NAME}_tran --rg-id $NORMAL_DATA_1_ID --rg PL:ILLUMINA --rg PU:${NORMAL_DATA_1_FC}-${NORMAL_DATA_1_BC}.${NORMAL_DATA_1_LN} --rg LB:${NORMAL_DATA_1_LB} --rg SM:${NORMAL_DATA_SM} --rna-strandness RF -1 $RNA_HOME/fastq/${NORMAL_DATA_1_BASE}_1.fastq -2 $RNA_HOME/fastq/${NORMAL_DATA_1_BASE}_2.fastq | $RNA_HOME/software/sambamba_v0.6.4 view -S -f bam -l 0 /dev/stdin | $RNA_HOME/software/sambamba_v0.6.4 sort -t $THREADS -m ${SAMBAMBA_RAM_GB}G --tmpdir $NORMAL_DATA_1_TEMP -o $RNA_HOME/alignments/${NORMAL_DATA_1_BASE}.bam /dev/stdin"

bsub -o $RNA_HOME/logs/hisat2_align_${NORMAL_DATA_2_BASE}.out -e $RNA_HOME/logs/hisat2_align_${NORMAL_DATA_2_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n $THREADS "$RNA_HOME/software/hisat2-2.0.4/hisat2 -p $THREADS --dta -x $RNA_HOME/refseq/${REFERENCE_NAME}_tran --rg-id $NORMAL_DATA_2_ID --rg PL:ILLUMINA --rg PU:${NORMAL_DATA_2_FC}-${NORMAL_DATA_2_BC}.${NORMAL_DATA_2_LN} --rg LB:${NORMAL_DATA_2_LB} --rg SM:${NORMAL_DATA_SM} --rna-strandness RF -1 $RNA_HOME/fastq/${NORMAL_DATA_2_BASE}_1.fastq -2 $RNA_HOME/fastq/${NORMAL_DATA_2_BASE}_2.fastq | $RNA_HOME/software/sambamba_v0.6.4 view -S -f bam -l 0 /dev/stdin | $RNA_HOME/software/sambamba_v0.6.4 sort -t $THREADS -m ${SAMBAMBA_RAM_GB}G --tmpdir $NORMAL_DATA_2_TEMP -o $RNA_HOME/alignments/${NORMAL_DATA_2_BASE}.bam /dev/stdin"

rmdir $TUMOR_DATA_1_TEMP/* $TUMOR_DATA_1_TEMP
rmdir $TUMOR_DATA_2_TEMP/* $TUMOR_DATA_2_TEMP

rmdir $NORMAL_DATA_1_TEMP/* $NORMAL_DATA_1_TEMP
rmdir $NORMAL_DATA_2_TEMP/* $NORMAL_DATA_2_TEMP

```

## Merge BAMs

```bash

# NORMAL

bsub -o $RNA_HOME/logs/sambamba_merge_${NORMAL_DATA_SM}.out -e $RNA_HOME/logs/sambamba_merge_${NORMAL_DATA_SM}.err -R "span[hosts=1]" -n $THREADS $RNA_HOME/software/sambamba_v0.6.4 merge -t $THREADS $RNA_HOME/alignments/${NORMAL_DATA_SM}.bam $RNA_HOME/alignments/${NORMAL_DATA_1_BASE}.bam $RNA_HOME/alignments/${NORMAL_DATA_2_BASE}.bam

#TUMOR

bsub -o $RNA_HOME/logs/sambamba_merge_${TUMOR_DATA_SM}.out -e $RNA_HOME/logs/sambamba_merge_${TUMOR_DATA_SM}.err -R "span[hosts=1]" -n $THREADS $RNA_HOME/software/sambamba_v0.6.4 merge -t $THREADS $RNA_HOME/alignments/${TUMOR_DATA_SM}.bam $RNA_HOME/alignments/${TUMOR_DATA_1_BASE}.bam $RNA_HOME/alignments/${TUMOR_DATA_2_BASE}.bam

```

## Transcript Assembly and Abundance Estimates

### Assemble Transcripts

```bash

bsub -o $RNA_HOME/logs/stringtie_assemble_${TUMOR_DATA_1_BASE}.out -e $RNA_HOME/logs/stringtie_assemble_${TUMOR_DATA_1_BASE}.err -M ${STRINGTIE_RAM_GB}000000 -R "select[mem>=${STRINGTIE_RAM_GB}000] rusage[mem=${STRINGTIE_RAM_GB}000] span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -G $REFERENCE_GTF -o $RNA_HOME/transcripts/${TUMOR_DATA_1_BASE}.gtf -p $THREADS -l $TUMOR_DATA_1_BASE $RNA_HOME/alignments/${TUMOR_DATA_1_BASE}.bam

bsub -o $RNA_HOME/logs/stringtie_assemble_${TUMOR_DATA_2_BASE}.out -e $RNA_HOME/logs/stringtie_assemble_${TUMOR_DATA_2_BASE}.err -M ${STRINGTIE_RAM_GB}000000 -R "select[mem>=${STRINGTIE_RAM_GB}000] rusage[mem=${STRINGTIE_RAM_GB}000] span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -G $REFERENCE_GTF -o $RNA_HOME/transcripts/${TUMOR_DATA_2_BASE}.gtf -p $THREADS -l $TUMOR_DATA_2_BASE $RNA_HOME/alignments/${TUMOR_DATA_2_BASE}.bam

bsub -o $RNA_HOME/logs/stringtie_assemble_${NORMAL_DATA_1_BASE}.out -e $RNA_HOME/logs/stringtie_assemble_${NORMAL_DATA_1_BASE}.err -M ${STRINGTIE_RAM_GB}000000 -R "select[mem>=${STRINGTIE_RAM_GB}000] rusage[mem=${STRINGTIE_RAM_GB}000] span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -G $REFERENCE_GTF -o $RNA_HOME/transcripts/${NORMAL_DATA_1_BASE}.gtf -p $THREADS -l $NORMAL_DATA_1_BASE $RNA_HOME/alignments/${NORMAL_DATA_1_BASE}.bam

bsub -o $RNA_HOME/logs/stringtie_assemble_${NORMAL_DATA_2_BASE}.out -e $RNA_HOME/logs/stringtie_assemble_${NORMAL_DATA_2_BASE}.err -M ${STRINGTIE_RAM_GB}000000 -R "select[mem>=${STRINGTIE_RAM_GB}000] rusage[mem=${STRINGTIE_RAM_GB}000] span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -G $REFERENCE_GTF -o $RNA_HOME/transcripts/${NORMAL_DATA_2_BASE}.gtf -p $THREADS -l $NORMAL_DATA_2_BASE $RNA_HOME/alignments/${NORMAL_DATA_2_BASE}.bam

```



#### Assemble Transcripts from Merged BAMs [OPTIONAL]

```base

bsub -o $RNA_HOME/logs/stringtie_assemble_${TUMOR_DATA_SM}.out -e $RNA_HOME/logs/stringtie_assemble_${TUMOR_DATA_SM}.err -M ${STRINGTIE_RAM_GB}000000 -R "select[mem>=${STRINGTIE_RAM_GB}000] rusage[mem=${STRINGTIE_RAM_GB}000] span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -G $REFERENCE_GTF -o $RNA_HOME/transcripts/${TUMOR_DATA_SM}.gtf -p $THREADS -l $TUMOR_DATA_SM $RNA_HOME/alignments/${TUMOR_DATA_SM}.bam

bsub -o $RNA_HOME/logs/stringtie_assemble_${NORMAL_DATA_SM}.out -e $RNA_HOME/logs/stringtie_assemble_${NORMAL_DATA_SM}.err -M ${STRINGTIE_RAM_GB}000000 -R "select[mem>=${STRINGTIE_RAM_GB}000] rusage[mem=${STRINGTIE_RAM_GB}000] span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -G $REFERENCE_GTF -o $RNA_HOME/transcripts/${NORMAL_DATA_SM}.gtf -p $THREADS -l $NORMAL_DATA_SM $RNA_HOME/alignments/${NORMAL_DATA_SM}.bam


```


### Merge Transcripts

```bash


$RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie --merge -p $LOCAL_THREADS -G $REFERENCE_GTF -o $RNA_HOME/transcripts/stringtie_merged.gtf $RNA_HOME/transcripts/${TUMOR_DATA_1_BASE}.gtf $RNA_HOME/transcripts/${TUMOR_DATA_2_BASE}.gtf $RNA_HOME/transcripts/${NORMAL_DATA_1_BASE}.gtf $RNA_HOME/transcripts/${NORMAL_DATA_2_BASE}.gtf

```

#### Merge Transcripts from Merged BAMs [OPTIONAL]

```base

$RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie --merge -p $LOCAL_THREADS -G $REFERENCE_GTF -o $RNA_HOME/transcripts/stringtie_merged_bams.gtf $RNA_HOME/transcripts/${TUMOR_DATA_SM}.gtf $RNA_HOME/transcripts/${NORMAL_DATA_SM}.gtf

```


### Compare Transcripts

```bash

$RNA_HOME/software/gffcompare-0.9.8.Linux_x86_64/gffcompare -r $REFERENCE_GTF -o $RNA_HOME/transcripts/gffcmp $RNA_HOME/transcripts/stringtie_merged.gtf

```

#### Compare Transcripts from Merged BAMs [OPTIONAL]

TODO: Add section to evaluate the transcripts before/after merging replicate BAMs.  Below should we simply compare to the merged GTF with replicates instead of $REFERENCE_GTF?
```bash

$RNA_HOME/software/gffcompare-0.9.8.Linux_x86_64/gffcompare -r $REFERENCE_GTF -o $RNA_HOME/transcripts/gffcmp2 $RNA_HOME/transcripts/stringtie_merged_bams.gtf

```

### Estimate Abundance

Continue with the merged GTF using all four replicates

```bash

mkdir -p $RNA_HOME/ballgown/${TUMOR_DATA_1_BASE}
mkdir -p $RNA_HOME/ballgown/${TUMOR_DATA_2_BASE}
mkdir -p $RNA_HOME/ballgown/${NORMAL_DATA_1_BASE}
mkdir -p $RNA_HOME/ballgown/${NORMAL_DATA_2_BASE}

bsub -o $RNA_HOME/logs/stringtie_abundance_${TUMOR_DATA_1_BASE}.out -e $RNA_HOME/logs/stringtie_abundance_${TUMOR_DATA_1_BASE}.err -R "span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -e -B -G $RNA_HOME/transcripts/gffcmp.annotated.gtf -o $RNA_HOME/ballgown/${TUMOR_DATA_1_BASE}/${TUMOR_DATA_1_BASE}.gtf -p $THREADS $RNA_HOME/alignments/${TUMOR_DATA_1_BASE}.bam

bsub -o $RNA_HOME/logs/stringtie_abundance_${TUMOR_DATA_2_BASE}.out -e $RNA_HOME/logs/stringtie_abundance_${TUMOR_DATA_2_BASE}.err -R "span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -e -B -G $RNA_HOME/transcripts/gffcmp.annotated.gtf -o $RNA_HOME/ballgown/${TUMOR_DATA_2_BASE}/${TUMOR_DATA_2_BASE}.gtf -p $THREADS $RNA_HOME/alignments/${TUMOR_DATA_2_BASE}.bam

bsub -o $RNA_HOME/logs/stringtie_abundance_${NORMAL_DATA_1_BASE}.out -e $RNA_HOME/logs/stringtie_abundance_${NORMAL_DATA_1_BASE}.err -R "span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -e -B -G $RNA_HOME/transcripts/gffcmp.annotated.gtf -o $RNA_HOME/ballgown/${NORMAL_DATA_1_BASE}/${NORMAL_DATA_1_BASE}.gtf -p $THREADS $RNA_HOME/alignments/${NORMAL_DATA_1_BASE}.bam

bsub -o $RNA_HOME/logs/stringtie_abundance_${NORMAL_DATA_2_BASE}.out -e $RNA_HOME/logs/stringtie_abundance_${NORMAL_DATA_2_BASE}.err -R "span[hosts=1]" -n $THREADS $RNA_HOME/software/stringtie-1.3.0.Linux_x86_64/stringtie -e -B -G $RNA_HOME/transcripts/gffcmp.annotated.gtf -o $RNA_HOME/ballgown/${NORMAL_DATA_2_BASE}/${NORMAL_DATA_2_BASE}.gtf -p $THREADS $RNA_HOME/alignments/${NORMAL_DATA_2_BASE}.bam

```


## Differential Expression

```bash

cd $RNA_HOME

# Make CSV/TSV for phenotype data
printf "\"ids\",\"type\"\n\"${NORMAL_DATA_1_BASE}\",\"normal\"\n\"${NORMAL_DATA_2_BASE}\",\"normal\"\n\"${TUMOR_DATA_1_BASE}\",\"tumor\"\n\"${TUMOR_DATA_2_BASE}\",\"tumor\"\n" > HCC1395.csv

$RNA_HOME/software/bin/R


```


```bash

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

# Load the phenotype data for each sample
pheno_data = read.csv("HCC1395.csv")

# Load ballgown data structures for each sample
bg = ballgown(dataDir = "ballgown", samplePattern = "H_NJ-HCC1395-HCC1395", pData=pheno_data)

# Filter low-abundance genes
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Identify signficant differently expressed Transcripts
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")

# Identify significant differently expressed Genes
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")

# Add gene names and gene IDs to the retuls_transcripts data frame
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)

# Sort from the smallest P value to largest
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

# Output as CSV
write.csv(results_transcripts,"HCC1395_transcript_results.csv",row.names=FALSE)
write.csv(results_genes,"HCC1395_genes_results.csv",row.names=FALSE)

# Output as TSV
write.table(results_transcripts,"HCC1395_transcript_results.tsv",sep="\t")
write.table(results_genes,"HCC1395_gene_results.tsv",sep="\t")

# Identify genes with p value < 0.05
subset(results_transcripts,results_transcripts$pval<0.05)
subset(results_genes,results_genes$pval<0.05)




```