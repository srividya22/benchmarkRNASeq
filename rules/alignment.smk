import os

def inputs(wildcards):
    """Returns fasta inputs for alignments."""
    in_path = reads_dir + "/sample_01" + {{pair}} + ".fasta"
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

def t_findex(wildcards,a_tool):
    """Returns the path to the transcriptome findex directory"""
    
    #tfindex_dir = expand("{base}/{org}/transcriptome/{aext}/{tool}_transcripts.idx".format(base=base_dir,org=org,aext=wildcards.aext,tool=a_tool))
    tfindex_dir = expand("{base}/{org}/transcriptome/{{aext}}/{tool}_transcripts.idx".format(base=base_dir,org=org,tool=a_tool))
    aexts = wildcards.aext
    return expand(gtf_path,tfindex_dir)

def get_gtf(wildcards):
    """Returns the path to the gtf file"""
 
    gtf_path = expand("{base}/{org}/annotations/gene_exons_{{aext}}.gtf".format(base=base_dir,org=org))
    aexts = wildcards.aext 
    return expand(gtf_path,aext=aexts)

def get_splicefile(wildcards):
    """Returns the path to the splice sites file"""
    
    splice_path = expand("{base}/{org}/annotations/splice_sites_{{aext}}.bed".format(base=base_dir,org=org))
    aexts=wildcards.aext
    return expand(splice_path,aext=aexts)

def get_logsdir(wildcards):
    """Return the path to the logs dir"""
   
    logs_path = expand("{lbase}/{{aext}}".format(lbase=logs_base))
    #aexts = list(config['A_EXTN'])
    aexts =  wildcards.aext
    return expand(logs_path,aext=aexts)
 
def get_outdir(wildcards):
    """ Returns the output directory to the analysis """
    outdir = expand("{rbase}/{{aext}}".format(rbase=results_base))
    aexts=  wildcards.aext
    return expand(outdir,aext=aexts)

print(inputs)
print(findex("bowtie2"))
print(get_gtf)
print(get_outdir)
print(get_logsdir)

rule tophat2:
     input: 
          #exe_dir=exe_dir,
          sample=inputs,
	  ref=findex("bowtie2"),
          gtf=get_gtf,          
          outdir=expand("{out}/tophat2".format(out=get_outdir))
     output:
          expand("{out}/tophat2/sample1/aligned_sorted.bam".format(out=get_outdir))
     log:
          logs_path=expand("{out}/tophat2".format(out=get_logsdir))
     resources:
          memory=10
     threads:
          config["tophat2"]["threads"]
     message: """--- running Tophat2"""
     shell:
            """ 
            echo "Sample : ${{input.sample}} , GTF : ${{input.gtf}} , Output : ${{input.outdir}}" 
            mkdir -p {input.outdir}
            exe_dir/tophat2 \
                   -o {input.outdir}/sample1 \
                   -G {input.gtf} \
                   -p {threads} \
                      {input.ref}/genome \
                      {input.sample} &> {input.outdir}/sample1/run.log
            cd {input.outdir}/sample1
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
          splice_bed=get_splicefile
          outdir=expand("{out}/hisat2_cuff".format(out=get_outdir))
     output:
          expand("{out}/hisat2_cuff/sample1/aligned_sorted.bam".format(out=get_outdir))
     log:
          logs_path=expand("{out}/hisat2_cuff".format(out=get_logsdir))
     resources:
          memory=10
     threads:
          config["hisat2"]["threads"]
     message: """--- running hisat2"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/hisat2 -p {threads} \
                   -S {input.outdir}/sample1/accepted_hits.sam \
                   --dta-cufflinks \
                   --known-splicesite-infile {input.splice_bed} \  
                   -x {input.ref}/genome \
                   -1 {input.sample} -2 {input.sample} &> {input.outdir}/sample1/run.log
            cd {input.outdir}/sample1
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
          splice_bed=get_splicefile
          outdir=expand("{out}/hisat2_str".format(out=get_outdir))
     output:
          expand("{out}/hisat2_str/sample1/aligned_sorted.bam".format(out=get_outdir))
     log:
          logs_path=expand("{out}/hisat2_str".format(out=get_logsdir))
     resources:
          memory=10
     threads:
          config["hisat2"]["threads"]
     message: """--- running hisat2"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/hisat2 -p {threads} \
                   -S {input.outdir}/sample1/accepted_hits.sam \
                   --dta \
                   --known-splicesite-infile {input.splice_bed} \
                   -x {input.ref}/genome \
                   -1 {input.sample} -2 {input.sample} &> {input.outdir}/sample1/run.log
            cd {input.outdir}/sample1
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
          gtf=get_gtf
          outdir=expand("{out}/star".format(out=get_outdir))
     output:
          expand("{out}/star/aligned_sorted.bam".format(out=get_outdir))
     log:
          logs_path=expand("{out}/star/".format(out=get_logsdir))
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
                   --outSAMtype BAM Unsorted --quantMode GeneCounts &> {input.outdir}/sample1/run.log
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
          gtf=get_gtf
          outdir=expand("{out}/star_2pass".format(out=get_outdir))
     output:
          expand("{out}/star_2pass/aligned_sorted.bam".format(out=get_outdir))
     log:
          logs_path=expand("{out}/star_2pass/".format(out=get_logsdir))
     resources:
          memory=10
     threads:
          config["star"]["threads"]
     message: """--- running STAR"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/STAR --runThreadN {threads} \
                   --twopassMode Basic \
                   --outFileNamePrefix {input.outdir}/sample1/sample1 \
                   --genomeDir {input.ref} \
                   --sjdbGTFfile {input.gtf} \
                   --readFilesIn {input.sample} \
                   --outFilterScoreMinOverLread 0.49 \
                   --outFilterMatchNminOverLread 0.49 --outSAMunmapped Within \
                   --outSAMprimaryFlag OneBestScore --outSAMstrandField intronMotif \
                   --outSAMtype BAM Unsorted --quantMode GeneCounts &> {input.outdir}/sample1/run.log
            cd {input.outdir}/sample1
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
          ref=t_findex("sailfish")
          outdir=expand("{out}/sailfish".format(out=get_outdir))
     output:
          expand("{out}/sailfish/quant.sf".format(out=get_outdir))
     log:
          logs_path=expand("{out}/sailfish/".format(out=get_logsdir))
     resources:
          memory=10
     threads:
          config["sailfish"]["threads"]
     message: """--- running Sailfish"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/sailfish quant --numBootstraps 100 \
                   -p {threads} \
                   -i {input.ref} \
                   -l "IU" \
                   -1 {input.sample} -2 {input.sample} \
                   -o {input.outdir}/sample1 &> {input.outdir}/sample1/run.log
           """

rule salmon:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=t_findex("salmon")
          outdir=expand("{out}/salmon".format(out=get_outdir))
     output:
          expand("{out}/sailfish/quant.sf".format(out=get_outdir))
     log:
          logs_path=expand("{out}/salmon/".format(out=get_logsdir))
     resources:
          memory=10
     threads:
          config["salmon"]["threads"]
     message: """--- running Salmon"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/sailfish quant --numBootstraps 100 \
                   -p {threads} \
                   -i {input.ref} \
                   -l "IU" \
                   -1 {input.sample} -2 {input.sample} \
                   -o {input.outdir}/sample1 &> {input.outdir}/sample1/run.log
           """

rule kallisto:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=t_findex("kallisto")
          outdir=expand("{out}/kallisto".format(out=get_outdir))
     output:
          expand("{out}/sailfish/abundance.tsv".format(out=get_outdir))
     log:
          logs_path=expand("{out}/kallisto/".format(out=get_logsdir))
     resources:
          memory=10
     threads:
          config["kallisto"]["threads"]
     message: """--- running Kallisto"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/kallisto quant -b 100 \
                   -t {threads} \
                   -i {input.ref} \
                   -l "IU" \
                   -o {input.outdir}/sample1 \
                   {input.sample} {input.sample} &> {input.outdir}/sample1/run.log
                   -o {input.outdir}/sample1
           """

#soapsplice -d <2BWT findex prefix> -1 <reads_a> -2 <reads_b> -r <length of reads> -I <insert size> -o <prefix of output files> [Other Options]

rule soapsplice:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=findex("soapsplice"),
     output:
          outdir=get_outdir("soapsplice")
     log:
          logs_path=get_logsdir("kallisto")
     resources:
          memory=10
     threads:
          config["soapsplice"]["threads"]
     message: """--- running soapsplice"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/soapsplice \
                   -d {input.ref} \
                   -1 {input.sample} -2 {input.sample} \
                   -r 100 \
                   -I 50 \
                   -o {input.outdir}/sample1 &> {input.outdir}/sample1/run.log
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
          outdir=get_outdir("mapsplice2")
     log:
          logs_path=get_logsdir("mapsplice2")
     resources:
          memory=10
     threads:
          config["mapsplice"]["threads"]
     message: """--- running mapsplice2"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/mapsplice.py -p {threads} \
                   -1 {input.sample} -2 {input.sample} \
                   -c {input.ref} \
                   -x {input.idx}\
                   -o {input.outdir}/sample1 &> {input.outdir}/sample1/run.log
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
          outdir = get_outdir("gsnap")
     log:
          logs_path=get_logsdir("gsnap")
     resources:
          memory=10
     threads:
          config[""]["threads"]
     message: """--- running mapsplice2"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/gsnap -t {threads} \
                   -D {input.ref}
                   -d {input.idx} \
                   {input.sample} {input.sample} \
                   -A sam > {input.outdir}/sample1/sample1.sam &> {input.outdir}/sample1/run.log
           """

#$ /PATH_TO_NOVOALIGN/novoalign -d /PATH_TO_PATHOGEN_GENOMES/genomes.nix -f /PATH_TO_SIMULATED_DATA/reads/simulated_set.fq -o SAM -r All > pathogen.sam

rule novoalign:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          idx=findex("genome.nix")
     output:
          outdir=get_outdir("novoalign")
     log:
          logs_path=get_logsdir("novoalign")
     resources:
          memory=10
     threads:
          config["novoalign"]["threads"]
     message: """--- running novoalign"""
     shell:"""
            mkdir -p {input.outdir}/sample1
            {input.exe_dir}/novoalign -p {threads} \
                   -d {input.idx} \
                   -f {input.sample} {input.sample} \
                   -o SAM -r All > {input.outdir}/sample1/sampl1.sam &> {input.outdir}/sample1/run.log
           """ 
