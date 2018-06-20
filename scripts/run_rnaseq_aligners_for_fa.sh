#!/bin/sh
## This Script run all the aligner jobs
exe_dir="$HOME/bin"

### Plant Species
org=$1

### Sample files
#sample_list=""sample_01" "sample_02" "sample_03" "sample_04""
sample_list=""sample_01" "sample_02""

### Set the Base directories
base_dir="/seq/schatz/sramakri/benchmark/RNASeq"
genome_base="${base_dir}/${org}/genome"
t_base="${base_dir}/${org}/transcriptome"
c_base="${base_dir}/${org}/contigs"
mt_base="${base_dir}/${org}/mutated_transcriptome"
gtf_base="${base_dir}/${org}/annotations"
input_base="${base_dir}/${org}/reads/"
results_base="${base_dir}/${org}/results"

### Experiment variables
### Run aligners for multiple coverages:  Low / Medium / High
covs=""5M" "10M" "15M" "20M" "25M""

### Run aligners mapping to multiple missing gene annotation set
parts=""0pcent" "5pcent" "10pcent" "20pcent" "25pcent""

### Run aligners mapping to  multiple heterozygosities
mutated=""0" "0.001" "0.005" "0.01" "0.02" "0.03" "0.05""

### Aligners Used in the Comparitive Study
aligners=""Tophat2" "Hisat2_stringtie" "Hisat2_cuff"  "STAR"  "STAR_2pass" "Sailfish" "Salmon" "kallisto" "kallisto-psuedo" "Salmon-quasi" "Sailfish-quasi" "Soapsplice" "GMAP" "Novoalign""


for i in ${parts}; do
	### Set the input variables
        #ref="${genome_base}/0/genome.fa"
        gtf="${gtf_base}/gene_exons_${i}.gtf"
        splice_file="${gtf_base}/splice_sites_${i}.bed"     
   
        ### Set of Output variables
        results_dir="${results_base}/alt_genome_annotations/${i}"
        logs_dir="${log_base}/alt_genome_annotations/${i}"
         
	
        for i in ${sample_list}; do
       
            file1="${input_base}/20M/${i}_1.fasta"
            file2="${input_base}/20M/${i}_2.fasta"
             
            ## Hisat2 Command
	    ref="${genome_base}/0/hisat2/genome"
            results_out=${results_dir}/Hisat2_stringtie
            mkdir -p ${results_out}
            
            cmd="$exe_dir/hisat2 -p $threads -S ${results_out}/Hisat2_stringtie/${i}.sam --dta --known-splicesite-infile $splice_file -x $ref -f -1 $file1 -2 $file2 && $HOME/bin/samtools view -bS -o ${results_out}/Hisat2_stringtie/${i}.bam ${results_out}/Hisat2_stringtie/${i}.sam && $HOME/bin/samtools view -h -F 4 -b -o ${results_out}/Hisat2_stringtie/${i}_allmapped.bam ${results_out}/Hisat2_stringtie/${i}.bam  && $HOME/bin/samtools sort -n -o ${results_out}/Hisat2_stringtie/${i}_nsorted.bam ${results_out}/Hisat2_stringtie/${i}_allmapped.bam && $HOME/bin/samtools sort -o ${results_out}/Hisat2_stringtie/${i}_sorted.bam ${results_out}/Hisat2_stringtie/${i}.bam && ${bed_dir}/bedtools bamtobed -i ${results_out}/Hisat2_stringtie/${i}_nsorted.bam > ${results_out}/Hisat2_stringtie/${i}.bed"

    	   echo ${cmd}
    	   jobid=`qsub -N align_hisat2_$i_$prefix ./jobscript.sh $cmd`
    	   echo "$i hisat2 $job_id" >> $results_out/aligner_jobids.out

	   ## Hisat2 to Cufflinks  Command

            results_out=${results_dir}/Hisat2_cuff
            mkdir -p ${results_out}
	    
            cmd="$exe_dir/hisat2 -p $threads -S ${results_out}/Hisat2_cuff/${i}.sam --dta-cufflinks --known-splicesite-infile $splice_file -x $ref -f -1 $file1 -2 $file2 && $HOME/bin/samtools view -bS -o ${results_out}/Hisat2_cuff/${i}.bam ${results_out}/Hisat2_cuff/${i}.sam && $HOME/bin/samtools view -h -F 4 -b -o ${results_out}/Hisat2_cuff/${i}_allmapped.bam ${results_out}/Hisat2_cuff/${i}.bam && $HOME/bin/samtools sort -n -o ${results_out}/Hisat2_cuff/${i}_nsorted.bam ${results_out}/Hisat2_cuff/${i}_allmapped.bam && $HOME/bin/samtools sort -o ${results_out}/Hisat2_cuff/${i}_sorted.bam ${results_out}/Hisat2_cuff/${i}.bam && ${bed_dir}/bedtools bamtobed -i ${results_out}//Hisat2_cuff/${i}_nsorted.bam > ${results_out}/Hisat2_cuff/${i}.bed"

    	    echo ${cmd}
            jobid=`qsub -N align_hisat2_$i_$prefix ./jobscript.sh $cmd`
            echo "$i hisat2 $jobid" >> $results_out/aligner_jobids.out

          ## Tophat Command
            
            ref="${genome_base}/0/bowtie2/genome"
            results_out=${results_dir}/Tophat2
            mkdir -p ${results_out}
 
	    cmd="${exe_dir}/tophat -p $threads -o ${results_out}/Tophat2/${i} -G $gtf $ref $file1 $file2 && $HOME/bin/samtools sort -n -o ${results_out}/Tophat2/${i}_nsorted.bam ${results_out}/Tophat2/${i}/accepted_hits.bam && $HOME/bin/samtools merge  ${results_out}/Tophat2/${i}_merged.bam ${results_out}/Tophat2/${i}/accepted_hits.bam ${results_out}/Tophat2/${i}/unmapped.bam && $HOME/bin/samtools sort -o ${results_out}/Tophat2/${i}_sorted.bam ${results_out}/Tophat2/${i}_merged.bam && ${bed_dir}/bedtools bamtobed -i ${results_out}/Tophat2/${i}_nsorted.bam > ${results_out}/Tophat2/${i}.bed"
	    
            echo ${cmd}
	    jobid=`qsub -N align_tophat2_$i_$prefix ./jobscript.sh $cmd`
	    echo "$i tophat2 $jobid" >> $base_out/aligner_jobids.out

	 ## STAR Command
           
            idx_dir="${genome_base}/0/"
            results_out=${results_dir}/STAR
            mkdir -p ${results_out}
	    
            cmd="${exe_dir}/STAR --runThreadN $threads --outFileNamePrefix ${results_out}/STAR/${i}_ --genomeDir $idx_dir --sjdbGTFfile $gtf --readFilesIn $file1 $file2 --outFilterScoreMinOverLread 0.49 --outFilterMatchNminOverLread 0.49 --outSAMunmapped Within --outSAMprimaryFlag OneBestScore --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --quantMode GeneCounts && $HOME/bin/samtools view -h -F 4 -b -o ${results_out}/STAR/${i}_allmapped.bam ${results_out}/STAR/${i}_Aligned.out.bam && $HOME/bin/samtools sort -n -o ${results_out}/STAR/${i}_nsorted.bam ${results_out}/STAR/${i}_allmapped.bam && $HOME/bin/samtools sort -o ${results_out}/STAR/${i}_sorted.bam ${results_out}/STAR/${i}_Aligned.out.bam && ${bed_dir}/bedtools bamtobed -i ${results_out}/STAR/${i}_allmapped.bam > ${results_out}/STAR/${i}.bed"
            
            echo ${cmd}
            jobid=`qsub -N align_star_$i_$prefix ./jobscript.sh $cmd`
            echo "$i STAR  $jobid" >> $results_out/aligner_jobids.out
	   
         ## STAR 2pass Command

            idx_dir="${genome_base}/0/"
            results_out=${results_dir}/STAR_2pass
            mkdir -p ${results_out}
	    
            cmd="${exe_dir}/STAR --runThreadN $threads --twopassMode Basic --outFileNamePrefix ${results_out}/STAR_2pass/${i}_ --genomeDir $idx_dir --sjdbGTFfile $gtf --readFilesIn $file1 $file2 --outFilterScoreMinOverLread 0.49 --outFilterMatchNminOverLread 0.49 --outSAMunmapped Within --outSAMprimaryFlag OneBestScore --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --quantMode GeneCounts && $HOME/bin/samtools view -h -F 4 -b -o ${results_out}/STAR_2pass/${i}_allmapped.bam ${results_out}/STAR_2pass/${i}_Aligned.out.bam  && $HOME/bin/samtools sort -n -o ${results_out}/STAR_2pass/${i}_nsorted.bam ${results_out}/STAR_2pass/${i}_allmapped.bam && $HOME/bin/samtools sort -o ${results_out}/STAR_2pass/${i}_sorted.bam ${results_out}/STAR_2pass/${i}_Aligned.out.bam && ${bed_dir}/bedtools bamtobed -i ${results_out}/STAR_2pass/${i}_allmapped.bam > ${results_out}/STAR_2pass/${i}.bed"

            echo ${cmd}
            jobid=`qsub -N align_star_$i_$prefix ./jobscript.sh $cmd`
            echo "$i STAR  $jobid" >> $results_out/aligner_jobids.out
 
        ## Sailfish Command
    
            ref="${genome_base}/0/"
            results_out=${results_dir}/STAR_2pass
            mkdir -p ${results_out}
            
            cmd="${exe_dir}/sailfish quant --numBootstraps 100 -p $threads -i $idx_dir -l "IU" -1 $file1 -2 $file2 -o ${base_out}/Sailfish/${i}"
    echo ${cmd}
             
            jobid=`qsub -N align_sailfish_$i_${prefix} ./jobscript.sh $cmd`
            echo "$i Sailfish $job_id" >> $base_out/aligner_jobids.out

    ## salmon  Command
    cmd="${HOME}/miniconda2/bin/salmon quant --numBootstraps 100 -p $threads -i ${idx_dir}/salmon_transcripts.idx -l "IU" -o ${base_out}/Salmon/${i} -1 $file1 -2 $file2"
    echo ${cmd}
    jobid=`qsub -N align_salmon_$i_${prefix} ./jobscript.sh $cmd`
    echo "$i Salmon $job_id" >> $base_out/aligner_jobids.out
#
#    ## Kallisto Command
#    cmd="${exe_dir}/kallisto quant -b 100 -t $threads -i ${idx_dir}/transcripts.idx -o ${base_out}/kallisto/${i} $file1 $file2"
#    echo ${cmd}
#    jobid=`qsub -N align_kallisto_$i ./jobscript.sh $cmd`
#    echo "$i Kallisto $job_id" >> $base_out/aligner_jobids.out

done                      
done

echo "##########"
echo "##########"
echo "Congratulations !! Finished Submitting Aligner Jobs.. Please check back for the results"
echo "##########"
echo "##########"
