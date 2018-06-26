#!/usr/bin/env bash
#submit_script.sh
                           # The following are options passed to qsub
#$ -V                      # Inherit the submission environment
#$ -cwd                    # Start job in submission directory
#$ -N snakemake-job
#$ -e $JOB_NAME.e$JOB_ID   # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID   # Name of the output file (eg. myMPI.oJobID)
#$ -m bes                  # Email at Begin and End of job or if suspended
#$ -l m_mem_free=10g
#$ -M srividya.ramki@gmail.com # E-mail address (change to your e-mail)

if [ "$#" -lt 1 ]; then
    echo "script.sh <out_dir (with config.yaml)> [additional snakemake arguments]*"
    exit
fi

source activate env-RNASeq
snakemake --use-conda -j 16 --local-cores 4 -w 90 \
    --cluster-config cluster.yml --cluster "qsub -k eo -m n -l nodes=1:ppn={cluster.n} -l mem={cluster.mem}gb -l walltime={cluster.time}" --directory "$@" 

# For Slurm jobs
#snakemake --use-conda --jobs 64 \
#    --cluster-config cluster.yml --cluster "sbatch --cpus-per-task={threads}" \
#    $@
