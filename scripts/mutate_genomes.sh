#!/bin/bash -e

#$ -V                      # Inherit the submission environment
#$ -cwd                    # Start job in submission directory
#$ -e $JOB_NAME.e$JOB_ID   # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID   # Name of the output file (eg. myMPI.oJobID)
#$ -m bes                  # Email at Begin and End of job or if suspended
#$ -pe mpi 8
#$ -l m_mem_free=10g
#$ -M srividya.ramki@gmail.com         # E-mail address (change to your e-mail)
#$ -N mutate_sequence        # The name of your job

exe="/seq/schatz/sramakri/benchmarkRNASeq/scripts/mutate_genomes.py"
ref=$1
outdir=$2
gtype=$3

prefix=$(basename $1 ".fa") 

mkdir -p ${outdir}/0
cp ${ref} ${outdir}/0/${gtype}.fa

for i in 0.001 0.005 0.01 0.02 0.03 0.05;
do
   mkdir -p ${outdir}/${i}
   cmd="python ${exe} ${ref} ${i} > ${outdir}/${i}/${gtype}.fa"
   qsub -N mutate_ref_${i} ./jobscript.sh ${cmd}
done
