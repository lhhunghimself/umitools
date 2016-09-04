#!/bin/bash
# This script converts a series of mRNA sequencing data file in FASTQ format
# to a table of UMI read counts of human genes in multiple sample conditions.


# 1 Parameters

# 1.1 Global

TOP_DIR="/mnt/backup/DetoxS/UMITest"

# 1.2 Dataset
SERIES="20150409"
SAMPLE_ID="RNAseq_${SERIES}"
LANES=6
DATA_DIR=$TOP_DIR
SEQ_DIR="${DATA_DIR}/Seqs"
ALIGN_DIR="${DATA_DIR}/output"
COUNT_DIR="${DATA_DIR}/Counts"
UMITOOLS_DIR="/mnt/backup/DetoxS/umitools/"
# 1.3 Reference
REF_DIR="$TOP_DIR/References/Broad_UMI"
SPECIES_DIR="${REF_DIR}/Human_RefSeq"
REF_SEQ_FILE="${SPECIES_DIR}/refMrna_ERCC_polyAstrip.hg19.fa"
SYM2REF_FILE="${SPECIES_DIR}/refGene.hg19.sym2ref.dat"
ERCC_SEQ_FILE="${REF_DIR}/ERCC92.fa"
BARCODE_FILE="${REF_DIR}/barcodes_trugrade_96_set4.dat"

# 1.4 Program
PROG_DIR="$TOP_DIR/Programs/Broad-DGE"
BWA_ALN_SEED_LENGTH=24
BWA_SAM_MAX_ALIGNS_FOR_XA_TAG=20
THREAD_NUMBER=8dif

# 2 Computation

# 2.1 Alignment
# Align sequence fragments to reference genome library.
let "IDX = 1"	
SEQ_FILES="";
#get files
while [ "$IDX" -le "${LANES}" ]; do
	SUBSAMPLE_ID="Lane$IDX"
	SEQ_FILE_R1="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R1.fastq"
	SEQ_FILE_R2="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R2.fastq"
	SEQ_FILES="${SEQ_FILES} ${SEQ_FILE_R1} ${SEQ_FILE_R2}"
	let "IDX = $IDX + 1"
done	
 #split into wells
 #use tight checking no mismatch no ambiguities to match original

echo "${UMITOOLS_DIR}/umisplit -v -l 16 -m 0 -N 0 -t $THREAD_NUMBER -b $BARCODE_FILE $SEQ_FILES"
$UMITOOLS_DIR/umisplit -v -l 16 -m 0 -N 0 -f -t $THREAD_NUMBER -b $BARCODE_FILE $SEQ_FILES

echo "${UMITOOLS_DIR}/scripts/multibwa.pl $THREAD_NUMBER >& /dev/null"
$UMITOOLS_DIR/scripts/multibwa.pl $THREAD_NUMBER >& /dev/null
echo "${UMITOOLS_DIR}/umimerge_parallel $SAMPLE_ID $SYM2REF_FILE $ERCC_SEQ_FILE $BARCODE_FILE $ALIGN_DIR $COUNT_DIR $THREAD_NUMBER"
$UMITOOLS_DIR/umimerge_parallel $SAMPLE_ID $SYM2REF_FILE $ERCC_SEQ_FILE $BARCODE_FILE $ALIGN_DIR $COUNT_DIR $THREAD_NUMBER
