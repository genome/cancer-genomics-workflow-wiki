
# Alignment

## Convert Unaligned BAMs to FASTQ

### Exome

Approximate Time : 10 Minutes
CPU Requested : 1
CPU Utlization : 80%
Max RAM : 512MB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/picard_sam-to-fastq_${TUMOR_EXOME_DATA_1_BASE}.out -e $SOMATIC_HOME/logs/picard_sam-to-fastq_${TUMOR_EXOME_DATA_1_BASE}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar SamToFastq I=$TUMOR_EXOME_DATA_1_BAM F=${SOMATIC_HOME}/fastq/${TUMOR_EXOME_DATA_1_BASE}_1.fastq F2=${SOMATIC_HOME}/fastq/${TUMOR_EXOME_DATA_1_BASE}_2.fastq

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/picard_sam-to-fastq_${TUMOR_EXOME_DATA_2_BASE}.out -e $SOMATIC_HOME/logs/picard_sam-to-fastq_${TUMOR_EXOME_DATA_2_BASE}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar SamToFastq I=$TUMOR_EXOME_DATA_2_BAM F=${SOMATIC_HOME}/fastq/${TUMOR_EXOME_DATA_2_BASE}_1.fastq F2=${SOMATIC_HOME}/fastq/${TUMOR_EXOME_DATA_2_BASE}_2.fastq

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/picard_sam-to-fastq_${NORMAL_EXOME_DATA_1_BASE}.out -e $SOMATIC_HOME/logs/picard_sam-to-fastq_${NORMAL_EXOME_DATA_1_BASE}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar SamToFastq I=$NORMAL_EXOME_DATA_1_BAM F=${SOMATIC_HOME}/fastq/${NORMAL_EXOME_DATA_1_BASE}_1.fastq F2=${SOMATIC_HOME}/fastq/${NORMAL_EXOME_DATA_1_BASE}_2.fastq

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/picard_sam-to-fastq_${NORMAL_EXOME_DATA_2_BASE}.out -e $SOMATIC_HOME/logs/picard_sam-to-fastq_${NORMAL_EXOME_DATA_2_BASE}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar SamToFastq I=$NORMAL_EXOME_DATA_2_BAM F=${SOMATIC_HOME}/fastq/${NORMAL_EXOME_DATA_2_BASE}_1.fastq F2=${SOMATIC_HOME}/fastq/${NORMAL_EXOME_DATA_2_BASE}_2.fastq

```

### WGS

Approximate Time :  Minutes
CPU Requested : 1
CPU Utlization : %
Max RAM : MB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/picard_sam-to-fastq_${TUMOR_WGS_DATA_1_BASE}.out -e $SOMATIC_HOME/logs/picard_sam-to-fastq_${TUMOR_WGS_DATA_1_BASE}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar SamToFastq I=$TUMOR_WGS_DATA_1_BAM F=${SOMATIC_HOME}/fastq/${TUMOR_WGS_DATA_1_BASE}_1.fastq F2=${SOMATIC_HOME}/fastq/${TUMOR_WGS_DATA_1_BASE}_2.fastq

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/picard_sam-to-fastq_${TUMOR_WGS_DATA_2_BASE}.out -e $SOMATIC_HOME/logs/picard_sam-to-fastq_${TUMOR_WGS_DATA_2_BASE}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar SamToFastq I=$TUMOR_WGS_DATA_2_BAM F=${SOMATIC_HOME}/fastq/${TUMOR_WGS_DATA_2_BASE}_1.fastq F2=${SOMATIC_HOME}/fastq/${TUMOR_WGS_DATA_2_BASE}_2.fastq

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/picard_sam-to-fastq_${NORMAL_WGS_DATA_1_BASE}.out -e $SOMATIC_HOME/logs/picard_sam-to-fastq_${NORMAL_WGS_DATA_1_BASE}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar SamToFastq I=$NORMAL_WGS_DATA_1_BAM F=${SOMATIC_HOME}/fastq/${NORMAL_WGS_DATA_1_BASE}_1.fastq F2=${SOMATIC_HOME}/fastq/${NORMAL_WGS_DATA_1_BASE}_2.fastq

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/picard_sam-to-fastq_${NORMAL_WGS_DATA_2_BASE}.out -e $SOMATIC_HOME/logs/picard_sam-to-fastq_${NORMAL_WGS_DATA_2_BASE}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar SamToFastq I=$NORMAL_WGS_DATA_2_BAM F=${SOMATIC_HOME}/fastq/${NORMAL_WGS_DATA_2_BASE}_1.fastq F2=${SOMATIC_HOME}/fastq/${NORMAL_WGS_DATA_2_BASE}_2.fastq

```

## BWA MEM

### Exome

Approximate Time : 25 Minutes
CPU Requested : 16
CPU Utlization : 80%
Max RAM : 6GB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/bwa_mem_samblaster_${TUMOR_EXOME_DATA_1_BASE}.out -e $SOMATIC_HOME/logs/bwa_mem_samblaster_${TUMOR_EXOME_DATA_1_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n ${THREADS} "bwa mem -K 10000000 -t ${THREADS} -Y -R \"@RG\tID:$TUMOR_EXOME_DATA_1_ID\tPL:ILLUMINA\tPU:${TUMOR_EXOME_DATA_1_FC}-${TUMOR_EXOME_DATA_1_BC}.${TUMOR_EXOME_DATA_1_LN}\tLB:${TUMOR_EXOME_DATA_1_LB}\tSM:${TUMOR_DATA_SM}\"  $SOMATIC_REFSEQ_FASTA ${SOMATIC_HOME}/fastq/${TUMOR_EXOME_DATA_1_BASE}_1.fastq ${SOMATIC_HOME}/fastq/${TUMOR_EXOME_DATA_1_BASE}_2.fastq | samblaster -a --addMateTags -o $SOMATIC_HOME/alignments/${TUMOR_EXOME_DATA_1_BASE}.sam"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/bwa_mem_samblaster_${TUMOR_EXOME_DATA_2_BASE}.out -e $SOMATIC_HOME/logs/bwa_mem_samblaster_${TUMOR_EXOME_DATA_2_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n ${THREADS} "bwa mem -K 10000000 -t ${THREADS} -Y -R \"@RG\tID:$TUMOR_EXOME_DATA_2_ID\tPL:ILLUMINA\tPU:${TUMOR_EXOME_DATA_2_FC}-${TUMOR_EXOME_DATA_2_BC}.${TUMOR_EXOME_DATA_2_LN}\tLB:${TUMOR_EXOME_DATA_2_LB}\tSM:${TUMOR_DATA_SM}\" $SOMATIC_REFSEQ_FASTA ${SOMATIC_HOME}/fastq/${TUMOR_EXOME_DATA_2_BASE}_1.fastq ${SOMATIC_HOME}/fastq/${TUMOR_EXOME_DATA_2_BASE}_2.fastq | samblaster -a --addMateTags -o $SOMATIC_HOME/alignments/${TUMOR_EXOME_DATA_2_BASE}.sam"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/bwa_mem_samblaster_${NORMAL_EXOME_DATA_1_BASE}.out -e $SOMATIC_HOME/logs/bwa_mem_samblaster_${NORMAL_EXOME_DATA_1_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n ${THREADS} "bwa mem -K 10000000 -t ${THREADS} -Y -R \"@RG\tID:$NORMAL_EXOME_DATA_1_ID\tPL:ILLUMINA\tPU:${NORMAL_EXOME_DATA_1_FC}-${NORMAL_EXOME_DATA_1_BC}.${NORMAL_EXOME_DATA_1_LN}\tLB:${NORMAL_EXOME_DATA_1_LB}\tSM:${NORMAL_DATA_SM}\" $SOMATIC_REFSEQ_FASTA ${SOMATIC_HOME}/fastq/${NORMAL_EXOME_DATA_1_BASE}_1.fastq ${SOMATIC_HOME}/fastq/${NORMAL_EXOME_DATA_1_BASE}_2.fastq | samblaster -a --addMateTags -o $SOMATIC_HOME/alignments/${NORMAL_EXOME_DATA_1_BASE}.sam"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/bwa_mem_samblaster_${NORMAL_EXOME_DATA_2_BASE}.out -e $SOMATIC_HOME/logs/bwa_mem_samblaster_${NORMAL_EXOME_DATA_2_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n ${THREADS} "bwa mem -K 10000000 -t ${THREADS} -Y -R \"@RG\tID:$NORMAL_EXOME_DATA_2_ID\tPL:ILLUMINA\tPU:${NORMAL_EXOME_DATA_2_FC}-${NORMAL_EXOME_DATA_2_BC}.${NORMAL_EXOME_DATA_2_LN}\tLB:${NORMAL_EXOME_DATA_2_LB}\tSM:${NORMAL_DATA_SM}\" $SOMATIC_REFSEQ_FASTA ${SOMATIC_HOME}/fastq/${NORMAL_EXOME_DATA_2_BASE}_1.fastq ${SOMATIC_HOME}/fastq/${NORMAL_EXOME_DATA_2_BASE}_2.fastq | samblaster -a --addMateTags -o $SOMATIC_HOME/alignments/${NORMAL_EXOME_DATA_2_BASE}.sam"

```

### WGS

Approximate Time :  24 Hours
CPU Requested : 8
CPU Utlization : 70%
Max RAM : 6 GB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/bwa_mem_samblaster_${TUMOR_WGS_DATA_1_BASE}.out -e $SOMATIC_HOME/logs/bwa_mem_samblaster_${TUMOR_WGS_DATA_1_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n ${THREADS} "bwa mem -K 10000000 -t ${THREADS} -Y -R \"@RG\tID:$TUMOR_WGS_DATA_1_ID\tPL:ILLUMINA\tPU:${TUMOR_WGS_DATA_1_FC}-${TUMOR_WGS_DATA_1_BC}.${TUMOR_WGS_DATA_1_LN}\tLB:${TUMOR_WGS_DATA_1_LB}\tSM:${TUMOR_DATA_SM}\"  $SOMATIC_REFSEQ_FASTA ${SOMATIC_HOME}/fastq/${TUMOR_WGS_DATA_1_BASE}_1.fastq ${SOMATIC_HOME}/fastq/${TUMOR_WGS_DATA_1_BASE}_2.fastq | samblaster -a --addMateTags -o $SOMATIC_HOME/alignments/${TUMOR_WGS_DATA_1_BASE}.sam"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/bwa_mem_samblaster_${TUMOR_WGS_DATA_2_BASE}.out -e $SOMATIC_HOME/logs/bwa_mem_samblaster_${TUMOR_WGS_DATA_2_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n ${THREADS} "bwa mem -K 10000000 -t ${THREADS} -Y -R \"@RG\tID:$TUMOR_WGS_DATA_2_ID\tPL:ILLUMINA\tPU:${TUMOR_WGS_DATA_2_FC}-${TUMOR_WGS_DATA_2_BC}.${TUMOR_WGS_DATA_2_LN}\tLB:${TUMOR_WGS_DATA_2_LB}\tSM:${TUMOR_DATA_SM}\" $SOMATIC_REFSEQ_FASTA ${SOMATIC_HOME}/fastq/${TUMOR_WGS_DATA_2_BASE}_1.fastq ${SOMATIC_HOME}/fastq/${TUMOR_WGS_DATA_2_BASE}_2.fastq | samblaster -a --addMateTags -o $SOMATIC_HOME/alignments/${TUMOR_WGS_DATA_2_BASE}.sam"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/bwa_mem_samblaster_${NORMAL_WGS_DATA_1_BASE}.out -e $SOMATIC_HOME/logs/bwa_mem_samblaster_${NORMAL_WGS_DATA_1_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n ${THREADS} "bwa mem -K 10000000 -t ${THREADS} -Y -R \"@RG\tID:$NORMAL_WGS_DATA_1_ID\tPL:ILLUMINA\tPU:${NORMAL_WGS_DATA_1_FC}-${NORMAL_WGS_DATA_1_BC}.${NORMAL_WGS_DATA_1_LN}\tLB:${NORMAL_WGS_DATA_1_LB}\tSM:${NORMAL_DATA_SM}\" $SOMATIC_REFSEQ_FASTA ${SOMATIC_HOME}/fastq/${NORMAL_WGS_DATA_1_BASE}_1.fastq ${SOMATIC_HOME}/fastq/${NORMAL_WGS_DATA_1_BASE}_2.fastq | samblaster -a --addMateTags -o $SOMATIC_HOME/alignments/${NORMAL_WGS_DATA_1_BASE}.sam"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/bwa_mem_samblaster_${NORMAL_WGS_DATA_2_BASE}.out -e $SOMATIC_HOME/logs/bwa_mem_samblaster_${NORMAL_WGS_DATA_2_BASE}.err -M ${ALIGN_RAM_GB}000000 -R "select[mem>=${ALIGN_RAM_GB}000] rusage[mem=${ALIGN_RAM_GB}000] span[hosts=1]" -n ${THREADS} "bwa mem -K 10000000 -t ${THREADS} -Y -R \"@RG\tID:$NORMAL_WGS_DATA_2_ID\tPL:ILLUMINA\tPU:${NORMAL_WGS_DATA_2_FC}-${NORMAL_WGS_DATA_2_BC}.${NORMAL_WGS_DATA_2_LN}\tLB:${NORMAL_WGS_DATA_2_LB}\tSM:${NORMAL_DATA_SM}\" $SOMATIC_REFSEQ_FASTA ${SOMATIC_HOME}/fastq/${NORMAL_WGS_DATA_2_BASE}_1.fastq ${SOMATIC_HOME}/fastq/${NORMAL_WGS_DATA_2_BASE}_2.fastq | samblaster -a --addMateTags -o $SOMATIC_HOME/alignments/${NORMAL_WGS_DATA_2_BASE}.sam"

```

## Merge

### Exome

Approximate Time : 25 Minutes
CPU Requested : 16
CPU Utlization : 10%
Max RAM : 256MB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_samtools_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_samtools_merge_${NORMAL_DATA_SM}.err -R "span[hosts=1]" -n $SORT_THREADS samtools merge -@ $SORT_THREADS $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_merged.bam $SOMATIC_HOME/alignments/${NORMAL_EXOME_DATA_1_BASE}.sam $SOMATIC_HOME/alignments/${NORMAL_EXOME_DATA_2_BASE}.sam

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_samtools_merge_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_samtools_merge_${TUMOR_DATA_SM}.err -R "span[hosts=1]" -n $SORT_THREADS samtools merge -@ $SORT_THREADS $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_merged.bam $SOMATIC_HOME/alignments/${TUMOR_EXOME_DATA_1_BASE}.sam $SOMATIC_HOME/alignments/${TUMOR_EXOME_DATA_2_BASE}.sam

```

### WGS

Approximate Time : 4 Hour
CPU Requested : 4
CPU Utlization : 50%
Max RAM : 70 MB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_samtools_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_samtools_merge_${NORMAL_DATA_SM}.err -R "span[hosts=1]" -n $SORT_THREADS samtools merge -@ $SORT_THREADS $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_merged.bam $SOMATIC_HOME/alignments/${NORMAL_WGS_DATA_1_BASE}.sam $SOMATIC_HOME/alignments/${NORMAL_WGS_DATA_2_BASE}.sam

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_samtools_merge_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_samtools_merge_${TUMOR_DATA_SM}.err -R "span[hosts=1]" -n $SORT_THREADS samtools merge -@ $SORT_THREADS $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_merged.bam $SOMATIC_HOME/alignments/${TUMOR_WGS_DATA_1_BASE}.sam $SOMATIC_HOME/alignments/${TUMOR_WGS_DATA_2_BASE}.sam

```

### Cleanup SAM files

```bash

find $SOMATIC_HOME/alignments/*.sam -exec echo rm {} \;

find $SOMATIC_HOME/alignments/*.sam -exec rm {} \;

```

## Queryname Sort

### Exome

Approximate Time : 20 Minutes
CPU Requested : 16
CPU Utlization : 40%
Max RAM : 16GB

```bash

TUMOR_EXOME_DATA_TEMP=`mktemp -d $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.XXXXXXXXXXXX`
NORMAL_EXOME_DATA_TEMP=`mktemp -d $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.XXXXXXXXXXXX`

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_sambamba_queryname_sort_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_sambamba_queryname_sort_${TUMOR_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n $SORT_THREADS sambamba sort -t $SORT_THREADS -m ${SAMBAMBA_RAM_GB}G -n --tmpdir=$TUMOR_EXOME_DATA_TEMP -o $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_queryname_sorted.bam $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_merged.bam

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_sambamba_queryname_sort_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_sambamba_queryname_sort_${NORMAL_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n $SORT_THREADS sambamba sort -t $SORT_THREADS -m ${SAMBAMBA_RAM_GB}G -n --tmpdir=$NORMAL_EXOME_DATA_TEMP -o $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_queryname_sorted.bam $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_merged.bam

rmdir $TUMOR_EXOME_DATA_TEMP/* $TUMOR_EXOME_DATA_TEMP
rmdir $NORMAL_EXOME_DATA_TEMP/* $NORMAL_EXOME_DATA_TEMP

rm $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_merged.bam
rm $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_merged.bam

```

### WGS

Approximate Time : Minutes
CPU Requested : 
CPU Utlization : %
Max RAM : GB

```bash

TUMOR_WGS_DATA_TEMP=`mktemp -d $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.XXXXXXXXXXXX`
NORMAL_WGS_DATA_TEMP=`mktemp -d $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.XXXXXXXXXXXX`

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_sambamba_queryname_sort_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_sambamba_queryname_sort_${TUMOR_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n $SORT_THREADS sambamba sort -t $SORT_THREADS -m ${SAMBAMBA_RAM_GB}G -n --tmpdir=$TUMOR_WGS_DATA_TEMP -o $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_queryname_sorted.bam $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_merged.bam

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_sambamba_queryname_sort_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_sambamba_queryname_sort_${NORMAL_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n $SORT_THREADS sambamba sort -t $SORT_THREADS -m ${SAMBAMBA_RAM_GB}G -n --tmpdir=$NORMAL_WGS_DATA_TEMP -o $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_queryname_sorted.bam $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_merged.bam

rmdir $TUMOR_WGS_DATA_TEMP/* $TUMOR_WGS_DATA_TEMP
rmdir $NORMAL_WGS_DATA_TEMP/* $NORMAL_WGS_DATA_TEMP

rm $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_merged.bam
rm $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_merged.bam

```

## Mark Duplicates and Position Sort

### Exome

Approximate Time : 45 Minutes
CPU Requested : 6
CPU Utlization : 40%
Max RAM : 28GB

```bash

TUMOR_EXOME_DATA_TEMP=`mktemp -d $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.XXXXXXXXXXXX`
NORMAL_EXOME_DATA_TEMP=`mktemp -d $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.XXXXXXXXXXXX`

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_picard_mrkdup_sambamba_sort_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_picard_mrkdup_sambamba_sort_${NORMAL_DATA_SM}.err -M ${MRKDUP_RAM_GB}000000 -R "select[mem>=${MRKDUP_RAM_GB}000] rusage[mem=${MRKDUP_RAM_GB}000] span[hosts=1]" -n ${MRKDUP_THREADS} " bash -c '$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar MarkDuplicates I=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_queryname_sorted.bam O=/dev/stdout ASSUME_SORT_ORDER=queryname METRICS_FILE=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | sambamba sort -t $SORT_THREADS -m ${SAMBAMBA_RAM_GB}G --tmpdir=$NORMAL_EXOME_DATA_TEMP -o $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_sorted.bam /dev/stdin' "

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_picard_mrkdup_sambamba_sort_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_picard_mrkdup_sambamba_sort_${TUMOR_DATA_SM}.err -M ${MRKDUP_RAM_GB}000000 -R "select[mem>=${MRKDUP_RAM_GB}000] rusage[mem=${MRKDUP_RAM_GB}000] span[hosts=1]" -n ${MRKDUP_THREADS} " bash -c '$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar MarkDuplicates I=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_queryname_sorted.bam O=/dev/stdout ASSUME_SORT_ORDER=queryname METRICS_FILE=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | sambamba sort -t $SORT_THREADS -m ${SAMBAMBA_RAM_GB}G --tmpdir=$TUMOR_EXOME_DATA_TEMP -o $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_sorted.bam /dev/stdin' "

rmdir $TUMOR_EXOME_DATA_TEMP/* $TUMOR_EXOME_DATA_TEMP
rmdir $NORMAL_EXOME_DATA_TEMP/* $NORMAL_EXOME_DATA_TEMP

rm $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_queryname_sorted.bam
rm $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_queryname_sorted.bam

```

### WGS

Approximate Time : Minutes
CPU Requested : 
CPU Utlization : %
Max RAM : GB

```bash

TUMOR_WGS_DATA_TEMP=`mktemp -d $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.XXXXXXXXXXXX`
NORMAL_WGS_DATA_TEMP=`mktemp -d $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.XXXXXXXXXXXX`

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_mrkdup_sambamba_sort_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_mrkdup_sambamba_sort_${NORMAL_DATA_SM}.err -M ${MRKDUP_RAM_GB}000000 -R "select[mem>=${MRKDUP_RAM_GB}000] rusage[mem=${MRKDUP_RAM_GB}000] span[hosts=1]" -n ${MRKDUP_THREADS} " bash -c '$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar MarkDuplicates I=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_queryname_sorted.bam O=/dev/stdout ASSUME_SORT_ORDER=queryname METRICS_FILE=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | sambamba sort -t $SORT_THREADS -m ${SAMBAMBA_RAM_GB}G --tmpdir=$NORMAL_WGS_DATA_TEMP -o $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_sorted.bam /dev/stdin' "

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_mrkdup_sambamba_sort_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_mrkdup_sambamba_sort_${TUMOR_DATA_SM}.err -M ${MRKDUP_RAM_GB}000000 -R "select[mem>=${MRKDUP_RAM_GB}000] rusage[mem=${MRKDUP_RAM_GB}000] span[hosts=1]" -n ${MRKDUP_THREADS} " bash -c '$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar MarkDuplicates I=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_queryname_sorted.bam O=/dev/stdout ASSUME_SORT_ORDER=queryname METRICS_FILE=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | sambamba sort -t $SORT_THREADS -m ${SAMBAMBA_RAM_GB}G --tmpdir=$TUMOR_WGS_DATA_TEMP -o $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_sorted.bam /dev/stdin' "

rmdir $TUMOR_WGS_DATA_TEMP/* $TUMOR_WGS_DATA_TEMP
rmdir $NORMAL_WGS_DATA_TEMP/* $NORMAL_WGS_DATA_TEMP

rm $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_queryname_sorted.bam
rm $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_queryname_sorted.bam

```

## Base Quality Score Recalibration (BQSR)

### Calculate BQSR Table

#### Exome

Approximate Time : 30 Minutes
CPU Requested : 4
CPU Utlization : 50%
Max RAM : 10GB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_gatk_base_recalibrator_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_gatk_base_recalibrator_${TUMOR_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n ${BQSR_THREADS} $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T BaseRecalibrator -R $SOMATIC_REFSEQ_FASTA -I $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_sorted.bam -o $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_bqsr.table -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_DBSNP}.gz -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_KNOWN_INDELS}  -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_MILLS_INDELS} --preserve_qscores_less_than 6 --disable_auto_index_creation_and_locking_when_reading_rods --disable_bam_indexing -dfrac .1 -nct $BQSR_THREADS -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_gatk_base_recalibrator_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_gatk_base_recalibrator_${NORMAL_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n ${BQSR_THREADS} $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T BaseRecalibrator -R $SOMATIC_REFSEQ_FASTA -I $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_sorted.bam -o $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_bqsr.table -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_DBSNP}.gz -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_KNOWN_INDELS}  -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_MILLS_INDELS} --preserve_qscores_less_than 6 --disable_auto_index_creation_and_locking_when_reading_rods --disable_bam_indexing -dfrac .1 -nct $BQSR_THREADS -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22

```
#### WGS

Approximate Time : Minutes
CPU Requested : 
CPU Utlization : %
Max RAM : GB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_gatk_base_recalibrator_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_gatk_base_recalibrator_${TUMOR_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n ${BQSR_THREADS} $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T BaseRecalibrator -R $SOMATIC_REFSEQ_FASTA -I $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_sorted.bam -o $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_bqsr.table -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_DBSNP}.gz -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_KNOWN_INDELS}  -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_MILLS_INDELS} --preserve_qscores_less_than 6 --disable_auto_index_creation_and_locking_when_reading_rods --disable_bam_indexing -dfrac .1 -nct $BQSR_THREADS -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_gatk_base_recalibrator_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_gatk_base_recalibrator_${NORMAL_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n ${BQSR_THREADS} $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T BaseRecalibrator -R $SOMATIC_REFSEQ_FASTA -I $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_sorted.bam -o $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_bqsr.table -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_DBSNP}.gz -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_KNOWN_INDELS}  -knownSites ${SOMATIC_REFSEQ_DIR}/${REFSEQ_MILLS_INDELS} --preserve_qscores_less_than 6 --disable_auto_index_creation_and_locking_when_reading_rods --disable_bam_indexing -dfrac .1 -nct $BQSR_THREADS -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22

```

### Apply BQSR

#### Exome

Approximate Time : 35 Minutes
CPU Requested : 8
CPU Utlization : 60%
Max RAM : 10GB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_gatk_print_reads_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_gatk_print_reads_${TUMOR_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n ${APPLY_BQSR_THREADS} $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T PrintReads -R $SOMATIC_REFSEQ_FASTA  -I $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_sorted.bam -o $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam -preserveQ 6 -BQSR $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_bqsr.table -SQQ 10 -SQQ 20 -SQQ 30 --disable_indel_quals -nct $APPLY_BQSR_THREADS

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_gatk_print_reads_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_gatk_print_reads_${NORMAL_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n ${APPLY_BQSR_THREADS} $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T PrintReads -R $SOMATIC_REFSEQ_FASTA  -I $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_sorted.bam -o $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.bam -preserveQ 6 -BQSR $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_bqsr.table -SQQ 10 -SQQ 20 -SQQ 30 --disable_indel_quals -nct $APPLY_BQSR_THREADS

```

#### WGS

Approximate Time :  Minutes
CPU Requested : 
CPU Utlization : %
Max RAM : GB

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_gatk_print_reads_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_gatk_print_reads_${TUMOR_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n ${APPLY_BQSR_THREADS} $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T PrintReads -R $SOMATIC_REFSEQ_FASTA  -I $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_sorted.bam -o $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam -preserveQ 6 -BQSR $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_bqsr.table -SQQ 10 -SQQ 20 -SQQ 30 --disable_indel_quals -nct $APPLY_BQSR_THREADS

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_gatk_print_reads_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_gatk_print_reads_${NORMAL_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000] span[hosts=1]" -n ${APPLY_BQSR_THREADS} $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T PrintReads -R $SOMATIC_REFSEQ_FASTA  -I $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_sorted.bam -o $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam -preserveQ 6 -BQSR $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_bqsr.table -SQQ 10 -SQQ 20 -SQQ 30 --disable_indel_quals -nct $APPLY_BQSR_THREADS

```

#### Cleanup Intermediate BAM files

```bash

find $SOMATIC_HOME/alignments/*_sorted.bam* -exec echo rm {} \;

find $SOMATIC_HOME/alignments/*_sorted.bam* -exec rm {} \;

```

# Quality Control

## SAMTools flagstat 

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_samtools_flagstat_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_samtools_flagstat_${NORMAL_DATA_SM}.err "samtools flagstat $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.bam > $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_flagstat.txt"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_samtools_flagstat_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_samtools_flagstat_${TUMOR_DATA_SM}.err "samtools flagstat $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam > $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_flagstat.txt"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_samtools_flagstat_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_samtools_flagstat_${NORMAL_DATA_SM}.err "samtools flagstat $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam > $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_flagstat.txt"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_samtools_flagstat_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_samtools_flagstat_${TUMOR_DATA_SM}.err "samtools flagstat $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam > $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_flagstat.txt"

```

## [Picard](https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics)

### Exome and WGS

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_picard_collect_insert_size_metrics_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_picard_collect_insert_size_metrics_${NORMAL_DATA_SM}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectInsertSizeMetrics I=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.bam O=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_insert_size_metrics.txt H=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_insert_size_metrics.pdf

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_picard_collect_insert_size_metrics_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_picard_collect_insert_size_metrics_${TUMOR_DATA_SM}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectInsertSizeMetrics I=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam O=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_insert_size_metrics.txt H=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_insert_size_metrics.pdf

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_collect_insert_size_metrics_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_collect_insert_size_metrics_${NORMAL_DATA_SM}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectInsertSizeMetrics I=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam O=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_insert_size_metrics.txt H=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_insert_size_metrics.pdf

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_collect_insert_size_metrics_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_collect_insert_size_metrics_${TUMOR_DATA_SM}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectInsertSizeMetrics I=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam O=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_insert_size_metrics.txt H=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_insert_size_metrics.pdf

```

## [Picard CollectAlignmentSummaryMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics)

### Exome and WGS

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_picard_collect_alignment_metrics_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_picard_collect_alignment_metrics_${NORMAL_DATA_SM}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectAlignmentSummaryMetrics I=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.bam O=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_alignment_metrics.txt R=$SOMATIC_REFSEQ_FASTA

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_picard_collect_alignment_metrics_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_picard_collect_alignment_metrics_${TUMOR_DATA_SM}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectAlignmentSummaryMetrics I=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam O=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_alignment_metrics.txt R=$SOMATIC_REFSEQ_FASTA

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_collect_alignment_metrics_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_collect_alignment_metrics_${NORMAL_DATA_SM}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectAlignmentSummaryMetrics I=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam O=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_alignment_metrics.txt R=$SOMATIC_REFSEQ_FASTA

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_collect_alignment_metrics_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_collect_alignment_metrics_${TUMOR_DATA_SM}.err $JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectAlignmentSummaryMetrics I=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam O=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_alignment_metrics.txt R=$SOMATIC_REFSEQ_FASTA

```

## [Picard CollectHsMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectHsMetrics)

### Exome ONLY

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_picard_collect_hs_metrics_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_picard_collect_hs_metrics_${NORMAL_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectHsMetrics I=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.bam O=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome_hs_metrics.txt R=$SOMATIC_REFSEQ_FASTA BI=$SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-probes.interval_list TI=$SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-targets.interval_list

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_picard_collect_hs_metrics_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/exome_picard_collect_hs_metrics_${TUMOR_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectHsMetrics I=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam O=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome_hs_metrics.txt R=$SOMATIC_REFSEQ_FASTA BI=$SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-probes.interval_list TI=$SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-targets.interval_list

```

## [Picard CollectGcBiasMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectGcBiasMetrics)

### WGS ONLY

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_collect_gc_bias_metrics_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_collect_gc_bias_metrics_${NORMAL_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectGcBiasMetrics I=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam O=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_gc_bias_metrics.txt R=$SOMATIC_REFSEQ_FASTA CHART=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_gc_bias_metrics.pdf S=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_gc_bias_summary.txt

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_collect_gc_bias_metrics_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_collect_gc_bias_metrics_${TUMOR_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectGcBiasMetrics I=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam O=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_gc_bias_metrics.txt R=$SOMATIC_REFSEQ_FASTA CHART=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_gc_bias_metrics.pdf S=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_gc_bias_summary.txt

```

## [Picard CollectWgsMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics)

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_collect_wgs_metrics_${NORMAL_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_collect_wgs_metrics_${NORMAL_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectWgsMetrics I=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam O=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs_metrics.txt R=$SOMATIC_REFSEQ_FASTA INTERVALS=$SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_autosomal.interval_list

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_picard_collect_wgs_metrics_${TUMOR_DATA_SM}.out -e $SOMATIC_HOME/logs/wgs_picard_collect_wgs_metrics_${TUMOR_DATA_SM}.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectWgsMetrics I=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam O=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs_metrics.txt R=$SOMATIC_REFSEQ_FASTA INTERVALS=$SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_autosomal.interval_list

```

## verifyBamID

# TODO: Run verifyBamId (limit to regions for exome)

