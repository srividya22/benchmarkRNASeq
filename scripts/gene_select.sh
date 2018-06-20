#!/bin/bash -e
##
# Script to select gene annotations from gtf 
#
# Author : Srividya Ramakrishnan
#
# Usage : "$0  <gtf_file> <output_dir>"
##

if [ "$#" -ne 3 ]; then
    echo "Usage : $0 <gtf_file> <output_dir>"
else
    exit 1
fi

gtf=$1
out_dir=$2
genes_file="${out_dir}/genes.txt"

## Make directories under out_dir
mkdir -p ${out_dir}
cp ${gtf} ${out_dir}/gene_exons_0pcent.gtf

## Select the genes from gff3 file
genes=($(ps -ef | awk '{print $NF}' ${gtf} | sort | uniq))

num_genes=${#genes[@]}

echo "Total num of genes : ${num_genes}"
for i in 5 10 15 20 25; do
    var1=$(echo "${i} * ${num_genes} /100" | bc)
    shuf -e -n ${var1} "${genes[@]}"  > ${out_dir}/excluded_${i}pcent_genes.txt
    awk 'NR==FNR{a[$1]; next; } !( $NF in a) {print $0}' ${out_dir}/excluded_${i}pcent_genes.txt ${gtf} > ${out_dir}/gene_exons_${i}pcent.gtf
done
