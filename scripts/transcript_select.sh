#!/bin/bash -e
##
# Script to select gene annotations from gtf 
#
# Author : Srividya Ramakrishnan
#
# Usage : "$0  <transcripts fasta> <output_dir>"
##

echo $#
if [ "$#" -ne 2 ]; then
    echo "Usage : $0 <transcripts fasta> <output_dir>"
    exit 1
fi

t_file=$1
out_dir=$2

## Make directories under out_dir
mkdir -p ${out_dir}/0pcent

cp ${t_file} ${out_dir}/0pcent/transcripts.fa

## Select the transcripts from fasta file
t_ids=($(ps -ef | grep ">" ${t_file} | sort | uniq))

num_ts=${#t_ids[@]}

echo "Total num of transcripts : ${num_ts}"
for i in 5 10 15 20 25; do
    mkdir -p ${out_dir}/${i}pcent
    var1=$(echo "${i} * ${num_ts} /100" | bc)
    shuf -e -n ${var1} "${t_ids[@]}" > ${out_dir}/excluded_${i}pcent_transcripts.txt
    $HOME/bin/filter_fasta.py ${out_dir}/excluded_${i}pcent_transcripts.txt ${t_file} ${out_dir}/${i}pcent/transcripts.fa
done
