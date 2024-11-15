#!/bin/bash
#PBS -N STAR_RSEM_416
#PBS -l nodes=1:ppn=20,walltime=48:00:00
#PBS -q batch
#PBS -d /projects/csna/rnaseq/Chesler_R01_DO_IVSA_Striatum/STAR_RSEM/
#PBS -t 1-416%24
#PBS -m e
#PBS -M michael.saul@jax.org

# Getting RSEM module version 1.3.0
module load rsem/1.3.0

# Running version 2.6.1 of STAR
STAR=/home/saulm/STAR/bin/Linux_x86_64/STAR

# Getting prefix information
PREFIX=`head -n $PBS_ARRAYID ../file_prefices.txt | tail -n 1`
echo "Using STAR to align files from "$PREFIX": `date`"

MOUSEGTF="/projects/csna/rnaseq/CCFounders_Sham_Cocaine/genome/Mus_musculus.GRCm38.94_ERCC.gtf"
OUTDIR="/projects/csna/rnaseq/Chesler_R01_DO_IVSA_Striatum/STAR_RSEM/"$PREFIX"/"
GENOMEDIR="/projects/csna/rnaseq/CCFounders_Sham_Cocaine/genome/"
echo ""
echo "Output Directory: "$OUTDIR
echo "Index Directory: "$GENOMEDIR
echo "Mouse Reference GTF: "$MOUSEGTF
echo ""

OUTBAM=$OUTDIR"Aligned.sortedByCoord.out.bam"
OUTTRANSBAM=$OUTDIR"Aligned.toTranscriptome.out.bam"
REFLOC="/projects/csna/rnaseq/CCFounders_Sham_Cocaine/genome/ref/Mus_musculus.GRCm38.dna_sm.primary_assembly_ERCC"
OUTQUANT=$OUTDIR"/"$PREFIX"_RSEM_quant"
mkdir $OUTDIR

FASTQ_1="/projects/csna/rnaseq/Chesler_R01_DO_IVSA_Striatum/LaneALL/"$PREFIX"_LaneALL_R1.fastq.gz"
FASTQ_2="/projects/csna/rnaseq/Chesler_R01_DO_IVSA_Striatum/LaneALL/"$PREFIX"_LaneALL_R2.fastq.gz"

echo "Forward Read FASTQ: "$FASTQ_1
echo "Reverse Read FASTQ: "$FASTQ_2
echo "Output Directory: "$OUTDIR
echo "Output BAM file: "$OUTBAM
echo "Output Transcriptome BAM file: "$OUTTRANSBAM
echo "Output Quantification Prefix: "$OUTQUANT
echo ""

# Running STAR
$STAR --runThreadN 20 \
--genomeDir $GENOMEDIR \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--quantTranscriptomeBAMcompression -1 \
--quantTranscriptomeBan IndelSoftclipSingleend \
--outFileNamePrefix $OUTDIR \
--outSAMattributes NH HI AS nM \
--readFilesIn $FASTQ_1 $FASTQ_2

echo "Using RSEM to EM quantify transcriptome of sample "$PREFIX": `date`"
echo ""

# Running RSEM
rsem-calculate-expression --bam --paired-end -p 20 \
--estimate-rspd --append-names \
$OUTTRANSBAM \
$REFLOC \
$OUTQUANT

rm $OUTTRANSBAM

echo "Alignment and transcriptome EM quantification of sample "$PREFIX" complete: `date`"
