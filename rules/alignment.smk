def inputs(wildcards):
    """Returns fasta inputs for alignments."""

    base_path = reads_dir + "/"+ sample_01{{pair}} + ".fasta"
    pairs = ["_1", "_2"] if is_paired else [""]
    #base_path = reads_dir + "/"+ sample_01{{pair}}.fasta".format(
    #    sample=wildcards.sample, lane=wildcards.lane)
    pairs = ["_1", "_2"] if is_paired else [""]

    return expand(base_path, pair=pairs)

def input(wildcards):
    """Returns fasta inputs for alignments."""

    base_path = reads_dir + "/"+ sample_01{{pair}} + ".fasta"
    pairs = ["_1", "_2"] if is_paired else [""]
    #base_path = reads_dir + "/"+ sample_01{{pair}}.fasta".format(
    #    sample=wildcards.sample, lane=wildcards.lane)
    pairs = ["_1", "_2"] if is_paired else [""]

    return expand(base_path, pair=pairs)

def get_contigs(wildcards):
    """Returns the path to the genome index directory."""

    contigs_dir = base_dir + "/" + org + "/contigs/0/"
    return contigs_dir

def index(a_tool):
    """Returns the path to the genome index directory."""
    index_dir = base_dir + "/" + org + "/genome/0/a_tool"
    return index_dir

def t_index(wildcards,a_tool):
    """Returns the path to the transcriptome index directory"""
    
    tindex_dir = expand("{base}/{org}/transcriptome/{aext}/{tool}_transcripts.idx".format(base=base_dir,org=org,aext=wildcards.aext,tool=a_tool))
    return tindex_dir

def get_gtf(wildcards):
    """Returns the path to the gtf file"""
 
    gtf_path = expand("{base}/{org}/annotations/gene_exons_{aext}.gtf".format(base=base_dir,org=org,aext=wildcards.aext))
    return gtf_path

def get_splicefile(wildcards):
    """Returns the path to the splice sites file"""
    
    splice_path = expand("{base}/{org}/annotations/splice_sites_{aext}.bed".format(base=base_dir,org=org,aext=wildcards.aext))
    return splice_path

def get_logsdir(wildcards,a_tool):
    """Return the path to the logs dir"""
   
    logs_path = expand("{lbase}/{aext}/{tool}".format(lbase=logs_base,aext=wildcards.aext,tool=a_tool))
    return logs_path
 
def get_outdir(wildcards,a_tool):
    """ Returns the output directory to the analysis """
    
    outdir = expand("{rbase}/{aext}/{tool}/".format(rbase=results_base,aext=wildcards.aext,tool=a_tool))

rule tophat2:
     input: 
          exe_dir=exe_dir,
          sample=inputs
     output:
          outdir=get_outdir("tophat2")
     log:
          logs_path=get_logsdir("tophat2")
     params:
	  ref = index("bowtie2"),
          gtf = get_gtf          
     resources:
          memory=10
     threads:
          config["tophat2"]["threads"]
     message: """--- running Tophat2"""
     shell:"""
            mkdir -p {output.outdir}
            {input.exe_dir}/tophat2 \
                   -o {output.outdir}/sample1 \
                   -G {params.gtf} \
                   -p {threads} \
                      {params.ref}/genome \
                      {input.sample} &> {output.outdir}/sample1/run.log
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
          sample=inputs
     output:
          outdir=get_outdir("hisat2_cuff")
     log:
          logs_path=get_logsdir("hisat2_cuff")
     params:
          ref=index("hisat2"),
          splice_bed=get_splicefile
     resources:
          memory=10
     threads:
          config["hisat2"]["threads"]
     message: """--- running hisat2"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/hisat2 -p {threads} \
                   -S {output.outdir}/sample1/accepted_hits.sam \
                   --dta-cufflinks \
                   --known-splicesite-infile {params.splice_bed} \  
                   -x {params.ref}/genome \
                   -1 {input.sample} -2 {input.sample} &> {output.outdir}/sample1/run.log
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
          sample=inputs
     output:
          outdir=get_outdir("hisat2_str")
     log:
          logs_path=get_logsdir("hisat2_str")
     params:
          ref=index("hisat2"),
          splice_bed=get_splicefile
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
                   --known-splicesite-infile {params.splice_bed} \
                   -x {params.ref}/genome \
                   -1 {input.sample} -2 {input.sample} &> {output.outdir}/sample1/run.log
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
          sample=inputs
     output:
          outdir=get_outdir("star")
     log:
          logs_path=get_logsdir("star")
     params:
          ref=index("star"),
          gtf=get_gtf
     resources:
          memory=10
     threads:
          config["star"]["threads"]
     message: """--- running STAR"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/STAR --runThreadN {threads} \
                   --outFileNamePrefix {output.outdir}/sample1/sample1 \
                   --genomeDir {params.ref} \
                   --sjdbGTFfile {params.gtf} \
                   --readFilesIn {input.sample} \
                   --outFilterScoreMinOverLread 0.49 			 \
                   --outFilterMatchNminOverLread 0.49 --outSAMunmapped Within \
                   --outSAMprimaryFlag OneBestScore --outSAMstrandField intronMotif \
                   --outSAMtype BAM Unsorted --quantMode GeneCounts &> {output.outdir}/sample1/run.log
            cd {output.outdir}/sample1
            mv sample1_Aligned.out.bam aligned.bam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam aligned_nsort.bed
           """

rule star_2pass:
     input:
          exe_dir=conda_dir,
          sample=inputs
     output:
          outdir=get_outdir("star_2pass")
     log:
          logs_path=get_logsdir("star_2pass")
     params:
          ref=index("star"),
          gtf=get_gtf
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
                   --genomeDir {params.ref} \
                   --sjdbGTFfile {params.gtf} \
                   --readFilesIn {input.sample} \
                   --outFilterScoreMinOverLread 0.49 \
                   --outFilterMatchNminOverLread 0.49 --outSAMunmapped Within \
                   --outSAMprimaryFlag OneBestScore --outSAMstrandField intronMotif \
                   --outSAMtype BAM Unsorted --quantMode GeneCounts &> {output.outdir}/sample1/run.log
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
          sample=inputs
     output:
          outdir=get_outdir("sailfish")
     log:
          logs_path=get_logsdir("sailfish")
     params:
          ref=t_index("sailfish")
     resources:
          memory=10
     threads:
          config["sailfish"]["threads"]
     message: """--- running Sailfish"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/sailfish quant --numBootstraps 100 \
                   -p {threads} \
                   -i {params.ref} \
                   -l "IU" \
                   -1 {input.sample} -2 {input.sample} \
                   -o {output.outdir}/sample1 &> {output.outdir}/sample1/run.log
           """

rule salmon:
     input:
          exe_dir=conda_dir,
          sample=inputs
     output:
          outdir=get_outdir("salmon")
     log:
          logs_path=get_logsdir("salmon")
     params:
          ref=t_index("salmon")
     resources:
          memory=10
     threads:
          config["salmon"]["threads"]
     message: """--- running Salmon"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/sailfish quant --numBootstraps 100 \
                   -p {threads} \
                   -i {params.ref} \
                   -l "IU" \
                   -1 {input.sample} -2 {input.sample} \
                   -o {output.outdir}/sample1 &> {output.outdir}/sample1/run.log
           """

rule kallisto:
     input:
          exe_dir=conda_dir,
          sample=inputs
     output:
          outdir=get_outdir("kallisto")
     log:
          logs_path=get_logsdir("kallisto")
     params:
          ref=t_index("kallisto")
     resources:
          memory=10
     threads:
          config["kallisto"]["threads"]
     message: """--- running Kallisto"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/kallisto quant -b 100 \
                   -t {threads} \
                   -i {params.ref} \
                   -l "IU" \
                   -o {output.outdir}/sample1 \
                   {input.sample} {input.sample} &> {output.outdir}/sample1/run.log
                   -o {output.outdir}/sample1
           """

#soapsplice -d <2BWT index prefix> -1 <reads_a> -2 <reads_b> -r <length of reads> -I <insert size> -o <prefix of output files> [Other Options]

rule soapsplice:
     input:
          exe_dir=conda_dir,
          sample=inputs
     output:
          outdir=get_outdir("soapsplice")
     log:
     params:
          ref=index("soapsplice"),
          logs_path=get_logsdir("soapsplice")
     resources:
          memory=10
     threads:
          config["soapsplice"]["threads"]
     message: """--- running soapsplice"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/soapsplice \
                   -d {params.ref} \
                   -1 {input.sample} -2 {input.sample} \
                   -r 100 \
                   -I 50 \
                   -o {output.outdir}/sample1 &> {output.outdir}/sample1/run.log
           """

#OUTPUT_DIR=/home/small_results/ 
#READ_FILE_END1=/home/small_1.fastq
#READ_FILE_END2=/home/small_2.fastq
#MAPSPLICE_DIR=/home/MapSplice-v2.2.0/
#REF_GENOME=/home/hg19_sequence/
#BOWTIE_INDEX=/home/hg19_index/humanchridx_M


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
          sample=inputs
     output:
          outdir=get_outdir("mapsplice2")
     log:
          logs_path=get_logsdir("mapsplice2")
     params:
          ref=get_contigs,
          idx=index("bowtie2")
     resources:
          memory=10
     threads:
          config["mapsplice"]["threads"]
     message: """--- running mapsplice2"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/mapsplice.py -p {threads} \
                   -1 {input.sample} -2 {input.sample} \
                   -c {params.ref} \
                   -x {param.idx}\
                   -o {output.outdir}/sample1 &> {output.outdir}/sample1/run.log
           """

#Paired : gsnap -d <genome> <fastq_file_1> <fastq_file_2> [<fastq_file_3> <fastq_file_4>...]
#Single  :gsnap -d <genome> --force-single-end <fastq_file_1> [<fastq_file_2>...]

rule gsnap:
     input:
          exe_dir = conda_dir,
          sample = inputs
     output:
          outdir = get_outdir("gsnap")
     log:
          logs_path=get_logsdir("gsnap")
     params:
          ref = get_contigs,
          idx = index("genome")
     resources:
          memory=10
     threads:
          config[""]["threads"]
     message: """--- running mapsplice2"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/gsnap -t {threads} \
                   -D {param.ref}
                   -d {params.idx} \
                   {input.sample} {input.sample} \
                   -A sam > {output.outdir}/sample1/sample1.sam &> {output.outdir}/sample1/run.log
           """

#$ /PATH_TO_NOVOALIGN/novoalign -d /PATH_TO_PATHOGEN_GENOMES/genomes.nix -f /PATH_TO_SIMULATED_DATA/reads/simulated_set.fq -o SAM -r All > pathogen.sam

rule novoalign:
     input:
          exe_dir=conda_dir,
          sample=inputs
     output:
          outdir=get_outdir("novoalign")
     log:
          logs_path=get_logsdir("novoalign")
     params:
          idx=index("genome.nix")
     resources:
          memory=10
     threads:
          config["novoalign"]["threads"]
     message: """--- running novoalign"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/novoalign -p {threads} \
                   -d {params.idx} \
                   -f {input.sample} {input.sample} \
                   -o SAM -r All > {output.outdir}/sample1/sampl1.sam &> {output.outdir}/sample1/run.log
           """ 
