import os

# Benchmarking Alignments 

def inputs(wildcards):
    """Returns fasta inputs for alignments."""
    in_path = "{rdir}/sample01{{pair}}.fasta".format(rdir=reads_dir)
    #in_path = reads_dir + "/sample01" + {{pair}} + ".fasta"
    pairs = ["_1", "_2"] if is_paired else [""]
    print(expand(in_path, pair=pairs))
    return expand(in_path, pair=pairs)

def input(wildcards):
    """Returns fasta inputs for alignments."""

    in_path = reads_dir + "/sample_01" + {{pairs}} + ".fasta"
    pairs = ["_1", "_2"] if is_paired else [""]
    #base_path = reads_dir + "/"+ sample_01{{pair}}.fasta".format(
    #    sample=wildcards.sample, lane=wildcards.lane)
    #pairs = ["_1", "_2"] if is_paired else [""]
    print(expand(in_path, pair=pairs))
    return expand(in_path, pair=pairs)

def get_contigs(wildcards):
    """Returns the path to the genome findex directory."""

    contigs_dir = base_dir + "/" + org + "/contigs/0/"
    print(contigs_dir)
    return contigs_dir

def findex(a_tool):
    """Returns the path to the genome findex directory."""
    findex_dir = base_dir + "/" + org + "/genome/0/" + a_tool
    print(findex_dir)
    return findex_dir

#def t_findex(wildcards,a_tool):
def t_findex(a_tool):
    """Returns the path to the transcriptome findex directory"""
    
    tfindex_dir = expand("{base}/{org}/transcriptome/{aext}/{tool}",base=base_dir,org=org,aext= config['A_EXTN'],tool=a_tool)
    #aexts = config['A_EXTN']
    #print(expand(tfindex_dir,aext=aexts))
    #aexts = wildcards.aext
    #return expand(tfindex_dir,aext=aexts)
    #return expand(gtf_path,tfindex_dir)
    return tfindex_dir

tfdir = "{base}/{org}/transcriptome".format(base=base_dir,org=org)
print(tfdir)

def get_gtf(wildcards):
    """Returns the path to the gtf file"""
 
    gtf_path = expand("{base}/{org}/annotations/gene_exons_{aext}.gtf",base=base_dir,org=org,aext=config['A_EXTN'])
    #aexts = wildcards.aext 
    #return expand(gtf_path,aext=aexts)
    return gtf_path

gtf_p = "{base}/{org}/annotations".format(base=base_dir,org=org)
print(gtf_p)

def get_splicefile(wildcards):
    """Returns the path to the splice sites file"""
    
    splice_path = expand("{base}/{org}/annotations/splice_sites_{aext}.bed",base=base_dir,org=org,aext=config['A_EXTN'])
    #aexts=wildcards.aext
    #return expand(splice_path,aext=aexts)
    return splice_path


def get_logsdir(wildcards):
    """Return the path to the logs dir"""
   
    logs_path = expand("{lbase}/{aext}",lbase=logs_base,aext=config['A_EXTN'])
    #aexts = config['A_EXTN']
    #aexts =  wildcards.aext
    #return expand(logs_path,aext=aexts)
    return logs_path

print(logs_base)
 
def get_outdir(wildcards):
    """ Returns the output directory to the analysis """
    outdir = expand("{rbase}/{aext}",rbase=results_base,aext=config['A_EXTN'])
    #aexts=  wildcards.aext
    #return expand(outdir,aext=aexts)
    return outdir

print(results_base)

print(inputs)
#print(findex("bowtie2"))
#print(get_gtf)
#print(get_outdir)
#print(get_logsdir)

rule tophat2:
     input: 
          exe_dir=exe_dir,
          sample=inputs,
	  ref=findex("bowtie2"),
          gtf=expand("{p_gtf}/gene_exons_{aext}.gtf",p_gtf=gtf_p,aext=config['A_EXTN'])          
     output:
          outdir = expand("{odir}/{aext}/tophat2",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/tophat2/sample1/run.log",odir=results_base,aext=config['A_EXTN']),
          out_bam = expand("{odir}/{aext}/tophat2/sample1/aligned_sorted.bam",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="tophat2.log")
     resources:
          memory=10
     threads:
          config["tophat2"]["threads"]
     message: """--- running Tophat2"""
     shell:
            """ 
            echo "Sample : ${{input.sample}} , GTF : ${{input.gtf}} , Output : ${{input.outdir}}" 
            mkdir -p {output.outdir}
            exe_dir/tophat2 \
                   -o {output.outdir}/sample1 \
                   -G {input.gtf} \
                   -p {threads} \
                      {input.ref}/genome \
                      {input.sample} &> {output.out_log}
            cd {output.outdir}/sample1
            mv accepted_hits.bam align.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools merge aligned.bam align.bam unmapped.bam 
            samtools sort -o aligned_sorted.bam aligned.bam 
            bedtools bamtobed -i aligned_nsort.bam aligned_nsort.bed
            """

rule hisat2_cuff:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=findex("hisat2"),
          splice_bed=expand("{p_gtf}/splice_sites_{aext}.bed",p_gtf=gtf_p,aext=config['A_EXTN'])          
     output:
          outdir = expand("{odir}/{aext}/hisat2_cuff",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/hisat2_cuff/sample1/run.log",odir=results_base,aext=config['A_EXTN']),
          out_bam = expand("{odir}/{aext}/hisat2_cuff/sample1/aligned_sorted.bam",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="hisat2_cuff.log")
     resources:
          memory=10
     threads:
          config["hisat2"]["threads"]
     message: """--- running hisat2"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            mkdir -p {logs.logs_path}
            {input.exe_dir}/hisat2 -p {threads} \
                   -S {output.outdir}/sample1/accepted_hits.sam \
                   --dta-cufflinks \
                   --known-splicesite-infile {input.splice_bed} \  
                   -x {input.ref}/genome \
                   -1 {input.sample} -2 {input.sample} &> {output.out_log}
            cd {output.outdir}/sample1
            samtools view -bS -o aligned.bam accepted_hits.sam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam aligned_nsort.bed
            """

rule hisat2_str:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=findex("hisat2"),
          splice_bed=expand("{p_gtf}/splice_sites_{aext}.bed",p_gtf=gtf_p,aext=config['A_EXTN'])          
     output:
          outdir = expand("{odir}/{aext}/hisat2_str",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/hisat2_str/sample1/run.log",odir=results_base,aext=config['A_EXTN']),
          out_bam = expand("{odir}/{aext}/hisat2_str/sample1/aligned_sorted.bam",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="hisat2_str.log")
     resources:
          memory=10
     threads:
          config["hisat2"]["threads"]
     message: """--- running hisat2"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/hisat2 -p {threads} \
                   -S {output.outdir}/sample1/accepted_hits.sam \
                   --dta \
                   --known-splicesite-infile {input.splice_bed} \
                   -x {input.ref}/genome \
                   -1 {input.sample} -2 {input.sample} &> {output.out_log}
            cd {output.outdir}/sample1
            samtools view -bS -o aligned.bam accepted_hits.sam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam aligned_nsort.bed
           """  
rule star:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=findex("star"),
          gtf=expand("{p_gtf}/gene_exons_{aext}.gtf",p_gtf=gtf_p,aext=config['A_EXTN'])
     output:
          outdir = expand("{odir}/{aext}/star",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/star/sample1/run.log",odir=results_base,aext=config['A_EXTN']),
          out_bam = expand("{odir}/{aext}/star/sample1/aligned_sorted.bam",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="star.log")
     resources:
          memory=10
     threads:
          config["star"]["threads"]
     message: """--- running STAR"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/STAR --runThreadN {threads} \
                   --outFileNamePrefix {input.outdir}/sample1/sample1 \
                   --genomeDir {input.ref} \
                   --sjdbGTFfile {input.gtf} \
                   --readFilesIn {input.sample} \
                   --outFilterScoreMinOverLread 0.49 			 \
                   --outFilterMatchNminOverLread 0.49 --outSAMunmapped Within \
                   --outSAMprimaryFlag OneBestScore --outSAMstrandField intronMotif \
                   --outSAMtype BAM Unsorted --quantMode GeneCounts &> {outdir.out_log}
            cd {input.outdir}/sample1
            mv sample1_Aligned.out.bam aligned.bam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam aligned_nsort.bed
           """

rule star_2pass:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=findex("star"),
          gtf=expand("{p_gtf}/gene_exons_{aext}.gtf",p_gtf=gtf_p,aext=config['A_EXTN'])
     output:
          outdir = expand("{odir}/{aext}/star_2pass",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/star_2pass/sample1/run.log",odir=results_base,aext=config['A_EXTN']),
          out_bam = expand("{odir}/{aext}/star_2pass/sample1/aligned_sorted.bam",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="star_2pass.log")
     resources:
          memory=10
     threads:
          config["star"]["threads"]
     message: """--- running STAR"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/STAR --runThreadN {threads} \
                   --twopassMode Basic \
                   --outFileNamePrefix {output.outdir}/sample1/sample1 \
                   --genomeDir {input.ref} \
                   --sjdbGTFfile {input.gtf} \
                   --readFilesIn {input.sample} \
                   --outFilterScoreMinOverLread 0.49 \
                   --outFilterMatchNminOverLread 0.49 --outSAMunmapped Within \
                   --outSAMprimaryFlag OneBestScore --outSAMstrandField intronMotif \
                   --outSAMtype BAM Unsorted --quantMode GeneCounts &> {output.out_log}
            cd {output.outdir}/sample1
            mv sample1_Aligned.out.bam aligned.bam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam aligned_nsort.bed
           """

rule sailfish:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=expand("{tdir}/{aext}/sailfish",tdir=tfdir,aext=config['A_EXTN'])
     output:
          outdir = expand("{odir}/{aext}/sailfish",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/sailfish/sample1/run.log",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="sailfish.log")
     resources:
          memory=10
     threads:
          config["sailfish"]["threads"]
     message: """--- running Sailfish"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/sailfish quant --numBootstraps 100 \
                   -p {threads} \
                   -i {input.ref} \
                   -l "IU" \
                   -1 {input.sample} -2 {input.sample} \
                   -o {output.outdir}/sample1 &> {output.out_log}
           """

rule salmon:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=expand("{tdir}/{aext}/salmon",tdir=tfdir,aext=config['A_EXTN'])
     output:
          outdir = expand("{odir}/{aext}/salmon",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/salmon/sample1/run.log",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="salmon.log")
     resources:
          memory=10
     threads:
          config["salmon"]["threads"]
     message: """--- running Salmon"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/sailfish quant --numBootstraps 100 \
                   -p {threads} \
                   -i {input.ref} \
                   -l "IU" \
                   -1 {input.sample} -2 {input.sample} \
                   -o {output.outdir}/sample1 &> {output.out_log}
           """

rule kallisto:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=expand("{tdir}/{aext}/kallisto/transcripts.idx",tdir=tfdir,aext=config['A_EXTN'])
     output:
          outdir = expand("{odir}/{aext}/kallisto",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/kallisto/sample1/run.log",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="kallisto.log")
     resources:
          memory=10
     threads:
          config["kallisto"]["threads"]
     message: """--- running Kallisto"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/kallisto quant -b 100 \
                   -t {threads} \
                   -i {input.ref} \
                   -l "IU" \
                   -o {output.outdir}/sample1 \
                   {input.sample} {input.sample} &> {output.out_log}
           """

#soapsplice -d <2BWT findex prefix> -1 <reads_a> -2 <reads_b> -r <length of reads> -I <insert size> -o <prefix of output files> [Other Options]

rule soapsplice:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=findex("soapsplice")
     output:
          outdir = expand("{odir}/{aext}/soapsplice",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/soapsplice/sample1/run.log",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="soapsplice.log")
     resources:
          memory=10
     threads:
          config["soapsplice"]["threads"]
     message: """--- running soapsplice"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/soapsplice \
                   -d {input.ref} \
                   -1 {input.sample} -2 {input.sample} \
                   -r 100 \
                   -I 50 \
                   -o {output.outdir}/sample1 &> {output.out_log}
           """

#OUTPUT_DIR=/home/small_results/ 
#READ_FILE_END1=/home/small_1.fastq
#READ_FILE_END2=/home/small_2.fastq
#MAPSPLICE_DIR=/home/MapSplice-v2.2.0/
#REF_GENOME=/home/hg19_sequence/
#BOWTIE_INDEX=/home/hg19_findex/humanchridx_M


#python $MAPSPLICE_DIR/mapsplice.py \
#       -1 $READ_FILE_END1 \
#       -2 $READ_FILE_END2 \
#       -c $REF_GENOME \
#       -x $BOWTIE_INDEX \
#       -p 8 \
#       -o $OUTPUT_DIR 2>log.txt

rule mapsplice2:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=get_contigs,
          idx=findex("bowtie2")
     output:
          outdir = expand("{odir}/{aext}/mapsplice2",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/mapsplice2/sample1/run.log",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="mapsplice2.log")
     resources:
          memory=10
     threads:
          config["mapsplice"]["threads"]
     message: """--- running mapsplice2"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/mapsplice.py -p {threads} \
                   -1 {input.sample} -2 {input.sample} \
                   -c {input.ref} \
                   -x {input.idx}\
                   -o {output.outdir}/sample1 &> {output.out_log}
           """

#Paired : gsnap -d <genome> <fastq_file_1> <fastq_file_2> [<fastq_file_3> <fastq_file_4>...]
#Single  :gsnap -d <genome> --force-single-end <fastq_file_1> [<fastq_file_2>...]

rule gsnap:
     input:
          exe_dir = conda_dir,
          sample = inputs,
          ref = get_contigs,
          idx = findex("genome")
     output:
          outdir = expand("{odir}/{aext}/gsnap",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/gsnap/sample1/run.log",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="gsnap.log")
     resources:
          memory=10
     threads:
          config["gmap"]["threads"]
     message: """--- running mapsplice2"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/gmap -t {threads} \
                   -D {input.ref}
                   -d {input.idx} \
                   {input.sample} {input.sample} \
                   -A sam > {output.outdir}/sample1/sample1.sam &> {output.out_log}
           """

#$ /PATH_TO_NOVOALIGN/novoalign -d /PATH_TO_PATHOGEN_GENOMES/genomes.nix -f /PATH_TO_SIMULATED_DATA/reads/simulated_set.fq -o SAM -r All > pathogen.sam

rule novoalign:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          idx=findex("genome.nix")
     output:
          outdir = expand("{odir}/{aext}/novoalign",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/novoalign/sample1/run.log",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/{tlog}",odir=results_base,aext=config['A_EXTN'],tlog="novoalign")
     resources:
          memory=10
     threads:
          config["novoalign"]["threads"]
     message: """--- running novoalign"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/novoalign -p {threads} \
                   -d {input.idx} \
                   -f {input.sample} {input.sample} \
                   -o SAM -r All > {output.outdir}/sample1/sampl1.sam &> {output.out_log}
           """

