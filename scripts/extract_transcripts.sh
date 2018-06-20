#!/bin/bash -e

######################################################
#
# Script to get the transcripts fasta from genome and gff3

# Usage : extract_transcripts.sh <genome fasta> <gff3  file> <outdir>
# Author : Srividya Ramakrishnan
# Affilation : Johns Hopkins University, MD
#
######################################################

if [ "$#" -ne 3 ]; then
    echo "Usage : extract_transcripts.sh <genome fasta> <gff3  file> <outdir>"
    exit 0
fi

genome=$1
gff3=$2
outdir=$3


awk '$3=="mRNA"{print $1"\t"$4"\t"$5"\t"$7"\t"$9}' ${gff3} | \
      awk -F";" '{print $1}' | sed 's/ID=//g' | \
      awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' > ${outdir}/transcripts.bed

bedtools getfasta -fi ${genome} -bed ${outdir}/transcripts.bed -s -name -fo ${outdir}/transcripts.fasta

echo "DONE"

### END ###
