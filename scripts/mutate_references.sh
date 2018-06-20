#!/bin/bash -e

#$ -V                      # Inherit the submission environment
#$ -cwd                    # Start job in submission directory
#$ -e $JOB_NAME.e$JOB_ID   # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID   # Name of the output file (eg. myMPI.oJobID)
#$ -m bes                  # Email at Begin and End of job or if suspended
#$ -pe mpi 8
#$ -l m_mem_free=10g
#$ -M srividya.ramki@gmail.com         # E-mail address (change to your e-mail)

exe="${HOME}/bin/msbar"
ref=$1
ref_base=$(basename $ref |  cut -d"." -f1)
output_dir=$2
base_count=$(${HOME}/bin/getlengths $ref | cut -f2 -d" " | paste -sd+ | bc)
p01=$(echo "0.001 * ${base_count}/1" | bc)
echo $p01
p05=$(echo "0.005 * ${base_count}/1" | bc)
echo $p05
p1=$(echo "0.01 * ${base_count}/1" | bc)
echo $p1
p2=$(echo "0.02 * ${base_count}/1" | bc)
echo $p2
p3=$(echo "0.03 * ${base_count}/1" | bc)
echo $p3
p5=$(echo "0.05 * ${base_count}/1" | bc)
echo $p5
cmd1="${exe} -sequence ${ref} -count ${p01} -point 4 -block 0 -codon 0 -outseq $output_dir/${ref_base}_0.1pcent_sub_mutations.fa"
cmd2="${exe} -sequence ${ref} -count ${p05} -point 4 -block 0 -codon 0 -outseq $output_dir/${ref_base}_0.5pcent_sub_mutations.fa"
cmd3="${exe} -sequence ${ref} -count ${p1} -point 4 -block 0 -codon 0 -outseq $output_dir/${ref_base}_1pcent_sub_mutations.fa"
cmd4="${exe} -sequence ${ref} -count ${p1} -point 4 -block 0 -codon 0 -outseq $output_dir/${ref_base}_2pcent_sub_mutations.fa"
cmd5="${exe} -sequence ${ref} -count ${p3} -point 4 -block 0 -codon 0 -outseq $output_dir/${ref_base}_3pcent_sub_mutations.fa"
cmd6="${exe} -sequence ${ref} -count ${p5} -point 4 -block 0 -codon 0 -outseq $output_dir/${ref_base}_5pcent_sub_mutations.fa"

### Submit all the cmds
jobid1=`qsub -N mutate_seq_0.1 ./jobscript.sh $cmd1| cut -f3 -d" "`
jobid2=`qsub -N mutate_seq_0.5 ./jobscript.sh $cmd2 | cut -f3 -d" "`
jobid3=`qsub -N mutate_seq_1 ./jobscript.sh $cmd3 | cut -f3 -d" "`
jobid4=`qsub -N mutate_seq_2 ./jobscript.sh $cmd4 | cut -f3 -d" "`
jobid5=`qsub -N mutate_seq_3 ./jobscript.sh $cmd5 | cut -f3 -d" "`
jobid6=`qsub -N mutate_seq_5 ./jobscript.sh $cmd6 | cut -f3 -d" "`
