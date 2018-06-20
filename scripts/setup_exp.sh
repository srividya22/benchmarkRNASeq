#!/bin/bash -e

###################################################
# Script to simulate RNASeq data for the experiment 
#
# Usage : setup_exp.sh <reference fasta> <reference gff3> <txt file containing samples>
#
# Author : Srividya Ramakrishnan
# Afflication : Johns Hopkins University
#

###################################################
SCRIPTS_PATH="/seq/schatz/sramakri/benchmarkRNASeq/scripts"

if [ "$#" -ne 6 ]; then
    echo "${0} <org> <reference fasta> <reference gtf> <transcripts fasta> <txt file containing fasta samples> <output experiment dir>"
    exit 0
fi

org=$1
ref=$2
gtf=$3
t_ref=$4
sample_text=$5
exp_dir=$6


cwd=$( pwd )
# Setup 1: Create required directories for the experiment
base_dir=${exp_dir}/${org}
mkdir -p ${base_dir}

# Make list of directories
cd ${base_dir}
mkdir -p annotations contigs genome logs mutated_transcriptome reads results seq transcriptome
cd ${cwd}

# Mutate genomes and create mutated reference files
out_dir=${base_dir}/genome
${SCRIPTS_PATH}/mutate_genomes.sh ${ref} ${out_dir} "genome"

# Mutate transcriptome and create mutated transcriptome files
out_dir=${base_dir}/mutated_transcriptome
${SCRIPTS_PATH}/mutate_genomes.sh ${t_ref} ${out_dir} "transcripts"

# Generate Missing annotation sets for gff3
out_dir=${base_dir}/annotations
${SCRIPTS_PATH}/gene_select.sh ${gtf} ${out_dir}

# Generate Missing transcripts from  fasta
out_dir=${base_dir}/transcriptome
${SCRIPTS_PATH}/transcript_select.sh ${t_ref} ${out_dir}
 
# Simulate RNASeq reads from real data
out_dir=${base_dir}/reads
echo ${out_dir}
