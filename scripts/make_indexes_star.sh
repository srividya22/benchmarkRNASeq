#!/bin/bash -e

#$ -V                      # Inherit the submission environment
#$ -cwd                    # Start job in submission directory
#$ -e $JOB_NAME.e$JOB_ID   # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID   # Name of the output file (eg. myMPI.oJobID)
#$ -m bes
#$ -l m_mem_free=25g
#$ -N generate_all_indexes      # The name of your job

### Set the input/output directories 
exe_dir="$HOME/bin"
base_dir="/seq/schatz/sramakri/benchmark/RNASeq"
org=$1
index_dir="${base_dir}/${org}/indexes"
genome_dir="${base_dir}/${org}/genome"
gtf_dir="${base_dir}/${org}/annotations"
t_dir="${base_dir}/${org}/transcriptome"
mt_dir="${base_dir}/${org}/mutated_transcriptome"
contigs_dir="${base_dir}/${org}/contigs"

### Set prefix for genome and transcripts
g_base="genome"
t_base="transcripts"

## File paths
gtf="${gtf_dir}/gene_exons_0pcent.gtf"
ref="${genome_dir}/0/genome.fa"
index_base="${index_dir}/0pcent/"
ref_fasta_dir="${contigs_dir}/0/*.fa"

## Experiment  : 1 Make indexes for original genomes with altered annotations

## Build Indexes for GTF independent alignment tools

#cmd="bowtie2-build --threads 8 ${ref} ${index_base}"
#echo $cmd
#job_id=$( qsub -N bowtie2-build_0 ./jobscript.sh ${cmd} )
#
#cmd="hisat2-build -p 8 $ref ${index_base}"
#echo $cmd
#job_id=$( qsub -N hisat2-build_0 ./jobscript.sh ${cmd} )
#
#cmd="gmap_build -d ${index_base} ${ref_fasta_dir}"
#echo $cmd
#job_id=$( qsub -N gmap_build_0 ./jobscript.sh ${cmd} )
#
##cmd="subread-buildindex -o ${index_base} ${ref_fasta_dir}"
##echo $cmd
##job_id=$( qsub ./jobscript.sh ${cmd} )
#
#cmd="2bwt-builder ${ref_fasta_dir} ${index_base}"
#echo $cmd
#job_id=$( qsub -N bwt-builder_0 ./jobscript.sh ${cmd} )
#
#cmd="novoindex ${index_base}.nix ${ref_fasta_dir}"
#echo $cmd
#job_id=$( qsub -N novoindex_0 ./jobscript.sh ${cmd} )

parts=""0pcent" "5pcent" "10pcent" "20pcent" "25pcent""
#parts=""0pcent""

for i in ${parts};
do
        ### Set input paths
        gtf="${gtf_dir}/gene_exons_${i}.gtf"
        
        STARTMP_DIR="${index_dir}/${i}/TMP"
        if [ -d "${STARTMP_DIR}" ]; then rm -rf ${STARTMP_DIR} ; fi
        cmd="STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ${index_dir}/${i}/ --genomeFastaFiles ${ref} --sjdbGTFfile $gtf --outTmpDir ${STARTMP_DIR}"
        echo $cmd
        job_id=$( qsub -N STARindex_0_${i} ./jobscript.sh ${cmd} )
	
#        cmd="${HOME}/bin/sailfish index -p 8 -t ${t_dir}/${i}/${t_base}.fa -o ${t_dir}/${i}/sailfish_transcripts.idx"
#        echo $cmd
#	job_id=$( qsub -N sailfishIndex_0_${i} ./jobscript.sh ${cmd} )
#     
#        cmd="salmon index -p 8 -t ${t_dir}/${i}/${t_base}.fa -i ${t_dir}/${i}/salmon_transcripts.idx"
#        echo $cmd
#        job_id=$( qsub -N salmonIndex_0_${i} ./jobscript.sh ${cmd} )
#        
#        cmd="kallisto index -i ${t_dir}/${i}/kallisto_transcripts.idx ${t_dir}/${i}/${t_base}.fa"
#        echo $cmd
#        job_id=$( qsub -N kallistoIndex_0_${i} ./jobscript.sh ${cmd} )
       
done

## Experiment : 2 Make indexes for the heterogyous genomes with same annotation

## Constants
gtf="${gtf_dir}/gene_exons_0pcent.gtf"
mutated=""0" "0.001" "0.005" "0.01" "0.02" "0.03" "0.05""
#mutated="0"

for i in ${mutated};
do
    ### Set variable input paths
        ref="${genome_dir}/${i}/genome.fa"
        ref_fasta_dir="${contigs_dir}/${i}/*.fa"
        index_base="${genome_dir}/${i}/"
        t_ref="${mt_dir}/${i}/transcripts.fa"

        ## Build bowtie2-build
#        mkdir -p "${index_base}/bowtie2/"
#        cmd="bowtie2-build --threads 8 ${ref} ${index_base}/bowtie2/genome"
#        echo $cmd
#        job_id=$( qsub -N bowtie2-build_${i} ./jobscript.sh ${cmd} )
#
        STARTMP_DIR="${index_base}/TMP"
        if [ -d "${STARTMP_DIR}" ]; then rm -rf ${STARTMP_DIR} ; fi
        cmd="STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ${index_base} --genomeFastaFiles ${ref} --sjdbGTFfile $gtf --outTmpDir ${STARTMP_DIR}"
        echo $cmd
        job_id=$( qsub -N STARindex_${i} ./jobscript.sh ${cmd} )

#        mkdir -p "${index_base}/hisat2/"
#        cmd="hisat2-build -p 8 $ref ${index_base}/hisat2/genome"
#        echo $cmd
#        job_id=$( qsub -N hisat2-build_${i} ./jobscript.sh ${cmd} )
#
#        cmd="${HOME}/bin/sailfish index -p 8 -t ${t_ref} -o ${mt_dir}/${i}/sailfish_transcripts.idx"
#        echo $cmd
#        job_id=$( qsub -N sailfishIndex_${i} ./jobscript.sh ${cmd} )
#
#        cmd="salmon index -p 8 -t ${t_ref} -i ${mt_dir}/${i}/salmon_transcripts.idx"
#        echo $cmd
#        job_id=$( qsub -N salmonIndex_${i} ./jobscript.sh ${cmd} )
#
#        cmd="kallisto index -i ${mt_dir}/${i}/kallisto_transcripts.idx ${t_ref}"
#        echo $cmd
#        job_id=$( qsub -N kallistoIndex_${i} ./jobscript.sh ${cmd} )
#
#        cmd="gmap_build -d ${index_base}/genome ${ref_fasta_dir}"
#        echo $cmd
#        job_id=$( qsub -N gmap_build_${i} ./jobscript.sh ${cmd} )
#
#        #cmd="subread-buildindex -o ${index_base} ${ref_fasta_dir}"
#        #echo $cmd
#        #job_id=$( qsub ./jobscript.sh ${cmd} )
#
#        cmd="2bwt-builder ${ref_fasta_dir} ${index_base}/genome"
#        echo $cmd
#        job_id=$( qsub -N bwt-builder_${i} ./jobscript.sh ${cmd} )
#        
#        cmd="novoindex ${index_base}/genome.nix ${ref_fasta_dir}"
#        echo $cmd
#        job_id=$( qsub -N novoindex ./jobscript.sh ${cmd} )
done
