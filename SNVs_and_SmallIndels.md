# Call Variants

## VarScan

* NOTE: Using a pileup per sample produced malformed VCF entries (ex. A/T).  Using a single pileup, `--mpileup 1`, for normal and tumor BAMS, the VCF entries look ok (ex. A,T).

### Exome

Approximate Time : 45 Minutes
CPU Requested : 2
CPU Utlization : 70%
Max RAM : 1GB

```bash

bsub -q $LSB_QUEUE -R 'span[hosts=1]' -n 2 -o $SOMATIC_HOME/logs/varscan_exome.out -e $SOMATIC_HOME/logs/varscan_exome.err bash -c "$JAVA_EIGHT -Xmx4g -jar $SOMATIC_HOME/software/VarScan.v2.4.2.jar somatic <($SOMATIC_HOME/software/bin/samtools mpileup -l $SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-targets.bed --no-BAQ -f $SOMATIC_REFSEQ_FASTA $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.bam $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam) $SOMATIC_HOME/varscan/exome --mpileup 1 --output-vcf"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/varscan_processSomatic_snp_exome.out -e $SOMATIC_HOME/logs/varscan_processSomatic_snp_exome.err $JAVA_EIGHT -Xmx4g -jar $SOMATIC_HOME/software/VarScan.v2.4.2.jar processSomatic $SOMATIC_HOME/varscan/exome.snp.vcf $SOMATIC_HOME/varscan/exome.snp
bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/varscan_processSomatic_indel_exome.out -e $SOMATIC_HOME/logs/varscan_processSomatic_indel_exome.err $JAVA_EIGHT -Xmx4g -jar $SOMATIC_HOME/software/VarScan.v2.4.2.jar processSomatic $SOMATIC_HOME/varscan/exome.indel.vcf $SOMATIC_HOME/varscan/exome.indel

find $SOMATIC_HOME/varscan -name *.vcf -exec bgzip -f {} \;
find $SOMATIC_HOME/varscan -name *.vcf.gz -exec tabix -f {} \;

$JAVA_EIGHT -Xmx4g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T VariantFiltration -R $SOMATIC_REFSEQ_FASTA --variant $SOMATIC_HOME/varscan/exome.snp.Somatic.vcf.gz --mask $SOMATIC_HOME/varscan/exome.snp.Somatic.hc.vcf.gz --maskName "processSomatic" --filterNotInMask -o $SOMATIC_HOME/varscan/exome.snp.Somatic.hc.filter.vcf.gz
$JAVA_EIGHT -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T VariantFiltration -R $SOMATIC_REFSEQ_FASTA --variant $SOMATIC_HOME/varscan/exome.indel.Somatic.vcf.gz --mask $SOMATIC_HOME/varscan/exome.indel.Somatic.hc.vcf.gz --maskName "processSomatic" --filterNotInMask -o $SOMATIC_HOME/varscan/exome.indel.Somatic.hc.filter.vcf.gz

bcftools concat -a -o $SOMATIC_HOME/varscan/exome.vcf.gz -O z $SOMATIC_HOME/varscan/exome.snp.Somatic.hc.filter.vcf.gz $SOMATIC_HOME/varscan/exome.indel.Somatic.hc.filter.vcf.gz

tabix -f $SOMATIC_HOME/varscan/exome.vcf.gz

```

### WGS

* Ran through Docker research-hpc queue

Approximate Time : 36 Hours
CPU Requested : 2
CPU Utlization : ?
Max RAM : ?

```bash

bsub -q $LSB_QUEUE -R 'span[hosts=1]' -n 2 -o $SOMATIC_HOME/logs/varscan_wgs.out -e $SOMATIC_HOME/logs/varscan_wgs.err bash -c "$JAVA_EIGHT -Xmx4g -jar $SOMATIC_HOME/software/VarScan.v2.4.2.jar somatic <($SOMATIC_HOME/software/bin/samtools mpileup --no-BAQ -f $SOMATIC_REFSEQ_FASTA $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam) $SOMATIC_HOME/varscan/wgs --mpileup 1 --output-vcf"

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/varscan_processSomatic_snp_wgs.out -e $SOMATIC_HOME/logs/varscan_processSomatic_snp_wgs.err $JAVA_EIGHT -Xmx4g -jar $SOMATIC_HOME/software/VarScan.v2.4.2.jar processSomatic $SOMATIC_HOME/varscan/wgs.snp.vcf $SOMATIC_HOME/varscan/wgs.snp
bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/varscan_processSomatic_indel_wgs.out -e $SOMATIC_HOME/logs/varscan_processSomatic_indel_wgs.err $JAVA_EIGHT -Xmx4g -jar $SOMATIC_HOME/software/VarScan.v2.4.2.jar processSomatic $SOMATIC_HOME/varscan/wgs.indel.vcf $SOMATIC_HOME/varscan/wgs.indel

find $SOMATIC_HOME/varscan -name *.vcf -exec bgzip -f {} \;
find $SOMATIC_HOME/varscan -name *.vcf.gz -exec tabix -f {} \;

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/varscan_processSomatic_snp_VariantFiltration_wgs.out -e $SOMATIC_HOME/logs/varscan_processSomatic_snp_VariantFiltration_wgs.err $JAVA_EIGHT -Xmx4g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T VariantFiltration -R $SOMATIC_REFSEQ_FASTA --variant $SOMATIC_HOME/varscan/wgs.snp.Somatic.vcf.gz --mask $SOMATIC_HOME/varscan/wgs.snp.Somatic.hc.vcf.gz --maskName "processSomatic" --filterNotInMask -o $SOMATIC_HOME/varscan/wgs.snp.Somatic.hc.filter.vcf.gz
bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/varscan_processSomatic_indel_VariantFiltration_wgs.out -e $SOMATIC_HOME/logs/varscan_processSomatic_indel_VariantFiltration_wgs.err $JAVA_EIGHT -Xmx4g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T VariantFiltration -R $SOMATIC_REFSEQ_FASTA --variant $SOMATIC_HOME/varscan/wgs.indel.Somatic.vcf.gz --mask $SOMATIC_HOME/varscan/wgs.indel.Somatic.hc.vcf.gz --maskName "processSomatic" --filterNotInMask -o $SOMATIC_HOME/varscan/wgs.indel.Somatic.hc.filter.vcf.gz

bcftools concat -a -o $SOMATIC_HOME/varscan/wgs.vcf.gz -O z $SOMATIC_HOME/varscan/wgs.snp.Somatic.hc.filter.vcf.gz $SOMATIC_HOME/varscan/wgs.indel.Somatic.hc.filter.vcf.gz

tabix -f $SOMATIC_HOME/varscan/wgs.vcf.gz

```

## Strelka

* Strelka is missing the 'GT' FORMAT tag which is required by GATK CombineVariants

### Exome

Approximate Time : 2 Hours
CPU Requested : 8
CPU Utlization : 28%
Max RAM : 4GB

```bash

$SOMATIC_HOME/software/strelka-2.7.1.centos5_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.bam --tumorBam=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam --referenceFasta=$SOMATIC_REFSEQ_FASTA --exome --runDir=$SOMATIC_HOME/strelka/exome

bsub -q $LSB_QUEUE -R 'span[hosts=1]' -n 8 -o $SOMATIC_HOME/logs/strelka_make_exome.out -e $SOMATIC_HOME/logs/strelka_make_exome.err $SOMATIC_HOME/strelka/exome/runWorkflow.py -m local -j 8

# WAIT FOR PRIOR STEP TO COMPLETE

zcat $SOMATIC_HOME/strelka/exome/results/variants/somatic.snvs.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > $SOMATIC_HOME/strelka/exome/results/variants/somatic.snvs.gt.vcf

zcat $SOMATIC_HOME/strelka/exome/results/variants/somatic.indels.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > $SOMATIC_HOME/strelka/exome/results/variants/somatic.indels.gt.vcf 

find $SOMATIC_HOME/strelka/exome/results/variants -name *.vcf -exec bgzip -f {} \;
find $SOMATIC_HOME/strelka/exome/results/variants -name *.vcf.gz -exec tabix -f {} \;
 
$SOMATIC_HOME/software/bin/bcftools concat -a -o $SOMATIC_HOME/strelka/exome.vcf.gz -O z $SOMATIC_HOME/strelka/exome/results/variants/somatic.snvs.gt.vcf.gz $SOMATIC_HOME/strelka/exome/results/variants/somatic.indels.gt.vcf.gz

tabix $SOMATIC_HOME/strelka/exome.vcf.gz

```

### WGS

* Ran in Docker queue research-hpc

Approximate Time : 4 Hours
CPU Requested : 8
CPU Utlization : ?
Max RAM : ?

```bash

$SOMATIC_HOME/software/strelka-2.7.1.centos5_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam=$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam --tumorBam=$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam --referenceFasta=$SOMATIC_REFSEQ_FASTA --runDir=$SOMATIC_HOME/strelka/wgs

bsub -q $LSB_QUEUE -R 'span[hosts=1]' -n 8 -o $SOMATIC_HOME/logs/strelka_make_wgs.out -e $SOMATIC_HOME/logs/strelka_make_wgs.err $SOMATIC_HOME/strelka/wgs/runWorkflow.py -m local -j 8

# WAIT FOR PRIOR STEP TO COMPLETE

zcat $SOMATIC_HOME/strelka/wgs/results/variants/somatic.snvs.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > $SOMATIC_HOME/strelka/wgs/results/variants/somatic.snvs.gt.vcf

zcat $SOMATIC_HOME/strelka/wgs/results/variants/somatic.indels.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > $SOMATIC_HOME/strelka/wgs/results/variants/somatic.indels.gt.vcf

find $SOMATIC_HOME/strelka/wgs/results/variants -name *.vcf -exec bgzip -f {} \;
find $SOMATIC_HOME/strelka/wgs/results/variants -name *.vcf.gz -exec tabix -f {} \;

$SOMATIC_HOME/software/bin/bcftools concat -a -o $SOMATIC_HOME/strelka/wgs.vcf.gz -O z $SOMATIC_HOME/strelka/wgs/results/variants/somatic.snvs.gt.vcf.gz $SOMATIC_HOME/strelka/wgs/results/variants/somatic.indels.gt.vcf.gz

tabix $SOMATIC_HOME/strelka/wgs.vcf.gz

```

## MuTect2

### Exome

Approximate Time : 4 Hours
CPU Requested : 1
CPU Utlization : 105%
Max RAM : 2GB

```bash

ls -1 $SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-targets/temp_*_of_$REFERENCE_SEQUENCE_CHUNKS/scattered.interval_list | xargs -n 1 dirname | xargs -n 1 basename | xargs -I CHUNK -n 1 bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/mutect_exome_CHUNK.out -e $SOMATIC_HOME/logs/mutect_exome_CHUNK.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T MuTect2 --disable_auto_index_creation_and_locking_when_reading_rods -R $SOMATIC_REFSEQ_FASTA -I:tumor $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam -I:normal $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.bam --dbsnp ${SOMATIC_REFSEQ_DIR}/${REFSEQ_DBSNP}.gz --cosmic $SOMATIC_REFSEQ_DIR/Cosmic_v79.dictsorted.vcf.gz -o $SOMATIC_HOME/mutect/exome/CHUNK.vcf.gz -L $SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-targets/CHUNK/scattered.interval_list

ls -1 $SOMATIC_HOME/refseq/GRCh38DH/xgen-exome-research-panel-targets/temp_*_of_$REFERENCE_SEQUENCE_CHUNKS/scattered.interval_list | xargs -n 1 dirname | xargs -n 1 basename | xargs -I CHUNK -n 1 echo $SOMATIC_HOME/mutect/exome/CHUNK.vcf.gz > $SOMATIC_HOME/mutect/exome/vcf.fof

$SOMATIC_HOME/software/bin/bcftools concat --allow-overlaps --remove-duplicates --file-list $SOMATIC_HOME/mutect/exome/vcf.fof --output-type z --output $SOMATIC_HOME/mutect/exome.vcf.gz

tabix $SOMATIC_HOME/mutect/exome.vcf.gz

```

### WGS

Approximate Time : 4 Hours
CPU Requested : 1
CPU Utlization : 105%
Max RAM : 6.5 GB

```bash

ls -1 $SOMATIC_REFSEQ_DIR/temp_*_of_$REFERENCE_SEQUENCE_CHUNKS/scattered.interval_list | xargs -n 1 dirname | xargs -n 1 basename | xargs -I CHUNK -n 1 bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/mutect_wgs_CHUNK.out -e $SOMATIC_HOME/logs/mutect_wgs_CHUNK.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T MuTect2 --disable_auto_index_creation_and_locking_when_reading_rods -R $SOMATIC_REFSEQ_FASTA -I:tumor $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam -I:normal $SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam --dbsnp ${SOMATIC_REFSEQ_DIR}/${REFSEQ_DBSNP}.gz --cosmic $SOMATIC_REFSEQ_DIR/Cosmic_v79.dictsorted.vcf.gz -o $SOMATIC_HOME/mutect/wgs/CHUNK.vcf.gz -L $SOMATIC_REFSEQ_DIR/CHUNK/scattered.interval_list

ls -1 $SOMATIC_HOME/refseq/GRCh38DH/scattered_calling_intervals/temp_*_of_$REFERENCE_SEQUENCE_CHUNKS/scattered.interval_list | xargs -n 1 dirname | xargs -n 1 basename | grep temp | xargs -I CHUNK -n 1 echo $SOMATIC_HOME/mutect/wgs/CHUNK.vcf.gz > $SOMATIC_HOME/mutect/wgs/vcf.fof

bcftools concat --allow-overlaps --remove-duplicates --file-list $SOMATIC_HOME/mutect/wgs/vcf.fof --output-type z --output $SOMATIC_HOME/mutect/wgs.vcf.gz

tabix $SOMATIC_HOME/mutect/wgs.vcf.gz

```

## Pindel

### Exome

```bash

printf "$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam\t400\tTUMOR\n$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_exome.bam\t400\tNORMAL\n" > $SOMATIC_HOME/pindel/exome.config

grep -v '^@' $SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-targets.interval_list | cut -f 1 | sort | uniq | xargs -n 1 basename | xargs -I CHR -n 1 bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_pindel_CHR.out -e $SOMATIC_HOME/logs/exome_pindel_CHR.err -M ${PINDEL_RAM_GB}000000 -R "select[mem>=${PINDEL_RAM_GB}000] rusage[mem=${PINDEL_RAM_GB}000] span[hosts=1]" -n 4 pindel -f $SOMATIC_REFSEQ_FASTA -i $SOMATIC_HOME/pindel/exome.config -c CHR -o $SOMATIC_HOME/pindel/exome_CHR -T 4 -w 20

cat $SOMATIC_HOME/pindel/exome_* | grep ChrID > $SOMATIC_HOME/pindel/exome.pindel

printf "input=$SOMATIC_HOME/pindel/exome.pindel\nvaf=0.1\ncov=20\nhom=6\npindel2vcf=$SOMATIC_HOME/software/bin/pindel2vcf\nreference=$SOMATIC_REFSEQ_FASTA\nreferencename=$GENOME_ASSEMBLY\nreferencedate=20161216\noutput=$SOMATIC_HOME/pindel/exome.raw.vcf\n" > $SOMATIC_HOME/pindel/filter_exome.config

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_pindel_somatic_filter.out -e $SOMATIC_HOME/logs/exome_pindel_somatic_filter.err /usr/bin/perl $SOMATIC_HOME/software/pindel-0.2.5b8/somatic_filter/somatic_indelfilter.pl $SOMATIC_HOME/pindel/filter_exome.config

bgzip -f $SOMATIC_HOME/pindel/exome.raw.vcf
tabix -f $SOMATIC_HOME/pindel/exome.raw.vcf.gz

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_pindel_select_variants.out -e $SOMATIC_HOME/logs/exome_pindel_select_variants.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T SelectVariants -R $SOMATIC_REFSEQ_FASTA --excludeFiltered --variant $SOMATIC_HOME/pindel/exome.raw.vcf.gz -L $SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-targets.interval_list -o $SOMATIC_HOME/pindel/exome.vcf.gz

tabix -f $SOMATIC_HOME/pindel/exome.vcf.gz

rm -fR $SOMATIC_HOME/pindel/exome_*

```

### WGS

```bash

printf "$SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam\t400\tTUMOR\n$SOMATIC_HOME/alignments/${NORMAL_DATA_SM}_wgs.bam\t400\tNORMAL\n" > $SOMATIC_HOME/pindel/wgs.config

grep -v '^@' $SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_calling_regions.interval_list | cut -f 1 | sort | uniq | xargs -I CHR -n 1 bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_pindel_CHR.out -e $SOMATIC_HOME/logs/wgs_pindel_CHR.err -M ${PINDEL_RAM_GB}000000 -R "select[mem>=${PINDEL_RAM_GB}000] rusage[mem=${PINDEL_RAM_GB}000] span[hosts=1]" -n 4 pindel -f $SOMATIC_REFSEQ_FASTA -i $SOMATIC_HOME/pindel/wgs.config -c CHR -o $SOMATIC_HOME/pindel/wgs_CHR -T 4 -w 20

cat $SOMATIC_HOME/pindel/wgs_* | grep ChrID > $SOMATIC_HOME/pindel/wgs.pindel

printf "input=$SOMATIC_HOME/pindel/wgs.pindel\nvaf=0.1\ncov=20\nhom=6\npindel2vcf=$SOMATIC_HOME/software/bin/pindel2vcf\nreference=$SOMATIC_REFSEQ_FASTA\nreferencename=$GENOME_ASSEMBLY\nreferencedate=20161216\noutput=$SOMATIC_HOME/pindel/wgs.vcf\n" > $SOMATIC_HOME/pindel/filter_wgs.config

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_pindel_somatic_filter.out -e $SOMATIC_HOME/logs/wgs_pindel_somatic_filter.err /usr/bin/perl $SOMATIC_HOME/software/pindel-0.2.5b8/somatic_filter/somatic_indelfilter.pl $SOMATIC_HOME/pindel/filter_wgs.config

bgzip -f $SOMATIC_HOME/pindel/wgs.vcf
tabix -f $SOMATIC_HOME/pindel/wgs.vcf.gz

rm -fR $SOMATIC_HOME/pindel/wgs_*

```

# Merge and Normalize Variants

## Combine Variants

* Merge ALL variants into a final VCF file that contains per caller info (like the GMS detailed.vcf)
* This is one of the "final" products of this workflow

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_gatk_combine_variants_unique.out -e $SOMATIC_HOME/logs/exome_gatk_combine_variants_unique.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T CombineVariants -R $SOMATIC_REFSEQ_FASTA -genotypeMergeOptions UNIQUIFY --variant:varscan $SOMATIC_HOME/varscan/exome.vcf.gz --variant:strelka $SOMATIC_HOME/strelka/exome.vcf.gz --variant:mutect $SOMATIC_HOME/mutect/exome.vcf.gz --variant:pindel $SOMATIC_HOME/pindel/exome.vcf.gz -o $SOMATIC_HOME/exome.unique.vcf.gz

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_gatk_combine_variants_merge.out -e $SOMATIC_HOME/logs/exome_gatk_combine_variants_merge.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T CombineVariants -R $SOMATIC_REFSEQ_FASTA -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan,strelka,pindel --variant:varscan $SOMATIC_HOME/varscan/exome.vcf.gz --variant:strelka $SOMATIC_HOME/strelka/exome.vcf.gz --variant:mutect $SOMATIC_HOME/mutect/exome.vcf.gz --variant:pindel $SOMATIC_HOME/pindel/exome.vcf.gz -o $SOMATIC_HOME/exome.merged.vcf

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_gatk_combine_variants_unique.out -e $SOMATIC_HOME/logs/wgs_gatk_combine_variants_unique.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T CombineVariants -R $SOMATIC_REFSEQ_FASTA -genotypeMergeOptions UNIQUIFY --variant:varscan $SOMATIC_HOME/varscan/wgs.vcf.gz --variant:strelka $SOMATIC_HOME/strelka/wgs.vcf.gz --variant:mutect $SOMATIC_HOME/mutect/wgs.vcf.gz --variant:pindel $SOMATIC_HOME/pindel/wgs.vcf.gz -o $SOMATIC_HOME/wgs.unique.vcf.gz

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_gatk_combine_variants_merge.out -e $SOMATIC_HOME/logs/wgs_gatk_combine_variants_merge.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T CombineVariants -R $SOMATIC_REFSEQ_FASTA -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan,strelka,pindel --variant:varscan $SOMATIC_HOME/varscan/wgs.vcf.gz --variant:strelka $SOMATIC_HOME/strelka/wgs.vcf.gz --variant:mutect $SOMATIC_HOME/mutect/wgs.vcf.gz --variant:pindel $SOMATIC_HOME/pindel/wgs.vcf.gz -o $SOMATIC_HOME/wgs.merged.vcf

```

## Normalize Variants

* TODO: Need to determine correct set of tools to normalize varaints


```bash

#bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_gatk_left_align.out -e $SOMATIC_HOME/logs/exome_gatk_left_align.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R $SOMATIC_REFSEQ_FASTA --variant $SOMATIC_HOME/exome.merged.vcf -o $SOMATIC_HOME/exome.merged.norm.vcf --splitMultiallelics

#bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_gatk_left_align.out -e $SOMATIC_HOME/logs/wgs_gatk_left_align.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R $SOMATIC_REFSEQ_FASTA --variant $SOMATIC_HOME/wgs.merged.vcf -o $SOMATIC_HOME/wgs.merged.norm.vcf --splitMultiallelics

```

# Filter Variants

## False-Positive Filter

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/fpfilter_exome.out -e $SOMATIC_HOME/logs/fpfilter_exome.err /usr/bin/perl $SOMATIC_HOME/software/fpfilter.pl --bam-file $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_exome.bam --sample TUMOR --reference $SOMATIC_REFSEQ_FASTA --bam-readcount /usr/bin/bam-readcount0.7 --vcf-file $SOMATIC_HOME/exome.merged.vcf --output $SOMATIC_HOME/exome.merged.fpfilter.vcf

# TODO: Handle MNPs before fpfilter

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/fpfilter_wgs.out -e $SOMATIC_HOME/logs/fpfilter_wgs.err /usr/bin/perl $SOMATIC_HOME/software/fpfilter.pl --bam-file $SOMATIC_HOME/alignments/${TUMOR_DATA_SM}_wgs.bam --sample TUMOR --reference $SOMATIC_REFSEQ_FASTA --bam-readcount /usr/bin/bam-readcount0.7 --vcf-file $SOMATIC_HOME/wgs.merged.vcf --output $SOMATIC_HOME/wgs.merged.fpfilter.vcf

bgzip -f $SOMATIC_HOME/exome.merged.vcf
bgzip -f $SOMATIC_HOME/exome.merged.fpfilter.vcf

bgzip -f $SOMATIC_HOME/wgs.merged.vcf
bgzip -f $SOMATIC_HOME/wgs.merged.fpfilter.vcf

tabix -f $SOMATIC_HOME/exome.merged.vcf.gz
tabix -f $SOMATIC_HOME/exome.merged.fpfilter.vcf.gz

tabix -f $SOMATIC_HOME/wgs.merged.vcf.gz
tabix -f $SOMATIC_HOME/wgs.merged.fpfilter.vcf.gz

```

## Hard Filter

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_gatk_select_variants.out -e $SOMATIC_HOME/logs/exome_gatk_select_variants.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T SelectVariants -R $SOMATIC_REFSEQ_FASTA --excludeFiltered --variant $SOMATIC_HOME/exome.merged.fpfilter.vcf.gz -o $SOMATIC_HOME/exome.merged.fpfilter.pass.vcf.gz

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_gatk_select_variants.out -e $SOMATIC_HOME/logs/wgs_gatk_select_variants.err -M ${SORT_RAM_GB}000000 -R "select[mem>=${SORT_RAM_GB}000] rusage[mem=${SORT_RAM_GB}000]" $JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/GenomeAnalysisTK.jar -T SelectVariants -R $SOMATIC_REFSEQ_FASTA --excludeFiltered --variant $SOMATIC_HOME/wgs.merged.fpfilter.vcf.gz -o $SOMATIC_HOME/wgs.merged.fpfilter.pass.vcf.gz

```

# Annotate Variants

* Only annotate PASS variants

```bash

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/exome_vep.out -e $SOMATIC_HOME/logs/exome_vep.err /usr/bin/perl $SOMATIC_HOME/software/ensembl-vep/vep.pl -i $SOMATIC_HOME/exome.merged.fpfilter.pass.vcf.gz --cache --dir $VEP_CACHE --format vcf --vcf --plugin Downstream --plugin Wildtype --symbol --terms SO --flag_pick -o $SOMATIC_HOME/exome.merged.fpfilter.pass.annotated.vcf.gz

bsub -q $LSB_QUEUE -o $SOMATIC_HOME/logs/wgs_vep.out -e $SOMATIC_HOME/logs/wgs_vep.err /usr/bin/perl $SOMATIC_HOME/software/ensembl-vep/vep.pl -i $SOMATIC_HOME/wgs.merged.fpfilter.pass.vcf.gz --cache --dir $VEP_CACHE --format vcf --vcf --plugin Downstream --plugin Wildtype --symbol --terms SO --flag_pick -o $SOMATIC_HOME/wgs.merged.fpfilter.pass.annotated.vcf.gz

```
