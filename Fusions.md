
# Initialize

## Environment

```bash

# define REFERENCE_NAME, REFERENCE_FASTA, REFERENCE_GTF, TUMOR_RNA_BAM, TUMOR_DNA_BAM, NORMAL_RNA_BAM

source /path/to/Fusions_config.sh 

FUSION_HOME=/gscuser/jwalker/git/HCC1395/integrate/${REFERENCE_NAME}
FUSION_RAM_GB=64

```

## Working Directories

```bash

mkdir -p $FUSION_HOME/software
mkdir -p $FUSION_HOME/logs
mkdir -p $FUSION_HOME/bwts
mkdir -p $FUSION_HOME/annotation

```

## Installs

### Integrate

```bash

cd $FUSION_HOME/software/

wget --no-check-certificate https://sourceforge.net/projects/integrate-fusion/files/INTEGRATE.0.2.6.tar.gz
tar -xzvf INTEGRATE.0.2.6.tar.gz
cd INTEGRATE_0_2_6/
mkdir INTEGRATE-build
cd INTEGRATE-build/
cmake ../Integrate/ -DCMAKE_BUILD_TYPE=release
make

$FUSION_HOME/software/INTEGRATE_0_2_6/INTEGRATE-build/bin/Integrate

```

### Integrate-NEO

```bash

cd $FUSION_HOME/software

git clone https://github.com/ChrisMaherLab/INTEGRATE-Neo.git
cd INTEGRATE-Neo/
cd INTEGRATE-Neo-V-1.0.0/
chmod +x install.sh
./install.sh -o $PWD

./integrate-neo.py

```

### gtfToGenePred

```bash

cd $FUSION_HOME/software

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
chmod +x gtfToGenePred

./gtfToGenePred

```

## Index

```bash

cd $FUSION_HOME

bsub -o $FUSION_HOME/logs/integrate_index.out -e $FUSION_HOME/logs/integrate_index.err -M 6000000 -R 'select[mem>=6000] rusage[mem=6000]' $FUSION_HOME/software/INTEGRATE_0_2_6/INTEGRATE-build/bin/Integrate mkbwt $REFERENCE_FASTA

```

## Annotation

```bash

cd $FUSION_HOME/annotation

$FUSION_HOME/software/gtfToGenePred -genePredExt -geneNameAsName2 $REFERENCE_GTF ${REFERENCE_NAME}.genePred
cut -f 1-10,12 ${REFERENCE_NAME}.genePred > ${REFERENCE_NAME}_tmp.txt

echo -e "#${REFERENCE_NAME}.ensGene.name\t${REFERENCE_NAME}.ensGene.chrom\t${REFERENCE_NAME}.ensGene.strand\t${REFERENCE_NAME}.ensGene.txStart\t${REFERENCE_NAME}.ensGene.txEnd\t${REFERENCE_NAME}.ensGene.cdsStart\t${REFERENCE_NAME}.ensGene.cdsEnd\t${REFERENCE_NAME}.ensGene.exonCount\t${REFERENCE_NAME}.ensGene.exonStarts\t${REFERENCE_NAME}.ensGene.exonEnds\t${REFERENCE_NAME}.ensemblToGeneName.value" > annot.enseml.${REFERENCE_NAME}.txt

cat ${REFERENCE_NAME}_tmp.txt >> annot.enseml.${REFERENCE_NAME}.txt

```

# Fusion Workflow

## Call Fusions

```bash

cd $FUSION_HOME

bsub -o $FUSION_HOME/logs/integrate_fusion.out -e $FUSION_HOME/logs/integrate_fusion.err -M ${FUSION_RAM_GB}000000 -R "select[mem>=${FUSION_RAM_GB}000] rusage[mem=${FUSION_RAM_GB}000]" $FUSION_HOME/software/INTEGRATE_0_2_6/INTEGRATE-build/bin/Integrate fusion $REFERENCE_FASTA $FUSION_HOME/annotation/annot.enseml.${REFERENCE_NAME}.txt $FUSION_HOME/bwts $TUMOR_RNA_BAM $TUMOR_RNA_BAM $TUMOR_DNA_BAM $NORMAL_RNA_BAM

```

## Annotate Fusions

```bash

cd $FUSION_HOME

$FUSION_HOME/software/INTEGRATE-Neo/INTEGRATE-Neo-V-1.0.0/fusionBedpeAnnotator --reference-file $REFERENCE_FASTA --gene-annotation-file $FUSION_HOME/annotation/annot.enseml.${REFERENCE_NAME}.txt --input-file $FUSION_HOME/fusions.bedpe --output-file $FUSION_HOME/fusions_annotated.bedpe

```