#!/bin/bash
source activate env-RNASeq
echo $PWD

snakemake --js $PWD/scripts/jobscript.sh\
    --printshellcmds\
    --restart-times 4\
    --cluster-config $PWD/cluster.yml\
    --jobname 'make.{jobid}.{rulename}'\
    --keep-going\
    --stats $PWD/snakemake.stats\
    --timestamp\
    --rerun-incomplete\
    -j 100\
    --cluster 'qsub -l walltime={cluster.time} -l mem={cluster.mem} -l vmem={cluster.mem} -l pmem={cluster.mem} -l nodes=1:ppn={cluster.cores} -o {cluster.logdir} -e {cluster.logdir}'
