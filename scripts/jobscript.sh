#!/bin/bash
#submit_script.sh
                           # The following are options passed to qsub
#$ -V                      # Inherit the submission environment
#$ -cwd                    # Start job in submission directory
#$ -e $JOB_NAME.e$JOB_ID   # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID   # Name of the output file (eg. myMPI.oJobID)
#$ -m bes                  # Email at Begin and End of job or if suspended
#$ -l m_mem_free=10g
#$ -M srividya.ramki@gmail.com         # E-mail address (change to your e-mail)
#source activate env-RNASeq
#source activate env-RNASeq
cmd=$*
SECONDS=0
#command block that takes time to complete...
#........
`$cmd`

#ENDTIME=$(date +%s)
