LSB_QUEUE='research-hpc'

JAVA_EIGHT=/gapp/x64linux/opt/java/jre/jre1.8.0_31/bin/java

INDEX_RAM_GB=8

ALIGN_RAM_GB=8
THREADS=16

SAMBAMBA_RAM_GB=18
SORT_RAM_GB=$((${SAMBAMBA_RAM_GB} + 2))

SORT_THREADS=2

JVM_RAM_GB=16

PICARD_THREADS=2

MRKDUP_RAM_GB=$((${SORT_RAM_GB} + ${JVM_RAM_GB}))
MRKDUP_THREADS=$((${SORT_THREADS} + ${PICARD_THREADS}))

BQSR_THREADS=4
APPLY_BQSR_THREADS=$((${BQSR_THREADS} * 2))

PINDEL_RAM_GB=30

GENOME_SPECIES=human
GENOME_ASSEMBLY=GRCh38DH
GENOME_BASENAME=GRCh38_full_analysis_set_plus_decoy_hla
GENOME_URI=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/$GENOME_BASENAME.fa

SOMATIC_HOME=/gscmnt/gc2764/cad/HCC1395/arvados

# See TUMOR_DATA_SM and NORMAL_DATA_SM in config.sh
#SOMATIC_NORMAL_SAMPLE=H_NJ-HCC1395-HCC1395_BL
#SOMATIC_TUMOR_SAMPLE=H_NJ-HCC1395-HCC1395

SOMATIC_REFSEQ_DIR=$SOMATIC_HOME/refseq/$GENOME_ASSEMBLY
SOMATIC_REFSEQ_FASTA=$SOMATIC_REFSEQ_DIR/$GENOME_BASENAME.fa
SOMATIC_REFSEQ_DICT=$SOMATIC_REFSEQ_DIR/$GENOME_BASENAME.dict

#GATK
BASE_REFSEQ=Homo_sapiens_assembly38

REFSEQ_DBSNP=$BASE_REFSEQ.dbsnp138.vcf
REFSEQ_KNOWN_INDELS=$BASE_REFSEQ.known_indels.vcf.gz
REFSEQ_MILLS_INDELS=Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

REFERENCE_SEQUENCE_CHUNKS=50

#/gscmnt/gc2142/techd/mutect/cosmic_v75/Cosmic.v75.dictsorted.vcf
#SOMATIC_COSMIC_VCF=

# VEP
VEP_CACHE=$SOMATIC_HOME/software/VEP_cache

export PATH="$SOMATIC_HOME/software/bin:$PATH"
export LSB_SUB_ADDITIONAL='docker(registry.gsc.wustl.edu/genome/genome_perl_environment)'