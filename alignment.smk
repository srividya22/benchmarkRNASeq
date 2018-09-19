import os
import subprocess, threading,traceback

#
# Benchmarking Alignments 

def inputs(wildcards):
    """Returns fasta inputs for alignments."""
    try:
        in_path = expand("{rdir}/sample01{{pair}}.fasta",rdir=reads_dir)
    #in_path = reads_dir + "/sample01" + {{pair}} + ".fasta"
        pairs = ["_1", "_2"] if is_paired else [""]
        print(expand(in_path, pair=pairs))
    except:
        print(traceback.format_exc())
        #raise Exception("Failed getting Input Fasta files")
    #if  not is_paired:
    return expand(in_path, pair=pairs)
    #else:
    #    fwd = expand(in_path,pair="_1")
    #    rev = expand(in_path,pair="_2")
    #    return (fwd , rev)


def input(wildcards):
    """Returns fasta inputs for alignments."""

    in_path = reads_dir + "/sample_01" + {{pairs}} + ".fasta"
    pairs = ["_1", "_2"] if is_paired else [""]
    #_path = reads_dir + "/"+ sample_01{{pair}}.fasta".format(
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
    try:
        findex_dir = base_dir + "/" + org + "/genome/0/" + a_tool
        print(findex_dir)
    except:
        print(traceback.format_exc())
        #raise Exception("Failed getting geneome Index for {0}".format(a_tool)) 
    return findex_dir

#def t_findex(wildcards,a_tool):
def t_findex(a_tool):
    """Returns the path to the transcriptome findex directory"""
    try:
       tfindex_dir = "{base}/{org}/transcriptome/{aext}/{tool}"
    #aexts = config['A_EXTN']
    #print(expand(tfindex_dir,aext=aexts))
    #aexts = wildcards.aext
    #return expand(tfindex_dir,aext=aexts)
    #return expand(gtf_path,tfindex_dir)
       print(tfindex_dir)
    except:
       print(traceback.format_exc())
       #raise Exception("Failed getting index for transcriptome {0}".format(a_tool))
    return tfindex_dir

tfdir = base_dir + "/" + org + "/transcriptome"
print(tfdir)

def get_gtf(wildcards):
    """Returns the path to the gtf file"""
    try:
       gtf_path = "{base}/{org}/annotations/gene_exons_{aext}.gtf"
    #gtf_path = expand("{}/{org}/annotations/gene_exons_{aext}.gtf",base=base_dir,org=org,aext=config['A_EXTN'])
    #aexts = wildcards.aext 
    #return expand(gtf_path,aext=aexts)
       print(gtf_path)
    except:
       print(traceback.format_exc())
       #raise Exception("Failed getting annotation gtf")
    return gtf_path

gtf_p = base_dir + "/" + org + "/annotations"
print(gtf_p)

def get_splicefile(wildcards):
    """Returns the path to the splice sites file"""
    try:
       splice_path = "{base}/{org}/annotations/splice_sites_{aext}.bed"
       print(splice_path)
    except:
       print(traceback.format_exc())
       #raise Exception("Failed getting splice file e")
    #aexts=wildcards.aext
    #return expand(splice_path,aext=aexts)
    return splice_path


def get_logsdir(wildcards):
    """Return the path to the logs dir"""
   
    logs_path = expand("{l}/{{aext}}",l=logs_base)
    print(logs_path)
    #aexts = config['A_EXTN']
    #aexts =  wildcards.aext
    #return expand(logs_path,aext=aexts)
    return logs_path

#print(logs_path)
 
def get_outdir(wildcards):
    """ Returns the output directory to the analysis """
    outdir = "{odir}/{aext}"
    #aexts=  wildcards.aext
    #return expand(outdir,aext=aexts)
    return outdir

print(findex("bowtie2"))
#print(get_gtf)
#print(get_outdir)
#print(get_logsdir)

rule tophat2:
     input: 
          exe_dir=exe_dir,
          #r1,r2=inputs,
          sample=inputs,
	  ref=findex("bowtie2"),
          #gtf=expand("{p_gtf}/gene_exons_{aext}.gtf",p_gtf=gtf_p,aext=config['A_EXTN'])          
          gtf=expand("{p_gtf}/gene_exons_{{aext}}.gtf",p_gtf=gtf_p)          
     output:
          outdir = "{odir}/{aext}/tophat2",
          out_log = "{odir}/{aext}/tophat2/sample1/run.log",
          out_bam = "{odir}/{aext}/tophat2/sample1/aligned_sorted.bam"
     log:
          logs_path= "{odir}/{aext}/tophat2.log"
     resources:
          memory=10
     threads:
          config["tophat2"]["threads"]
     message: """--- running Tophat2"""
     shell:
            """
            source deactivate
            echo "Sample : ${input.sample} , GTF : ${input.gtf} , Output : ${output.outdir}" 
            mkdir -p {output.outdir}
            {input.exe_dir}/tophat2 \
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
            bedtools bamtobed -i aligned_nsort.bam > aligned_nsort.bed
            cd -
            source activate env-RNASeq
            """

rule hisat2_cuff:
     input:
          exe_dir=conda_dir,
          #r1,r2=inputs,
          sample=inputs,
          ref=findex("hisat2"),
          #splice_bed=expand("{p_gtf}/splice_sites_{aext}.bed",p_gtf=gtf_p,aext=config['A_EXTN'])          
          splice_bed=expand("{p_gtf}/splice_sites_{{aext}}.bed",p_gtf=gtf_p)          
     output:
          outdir = "{odir}/{aext}/hisat2_cuff",
          out_log = "{odir}/{aext}/hisat2_cuff/sample1/run.log",
          out_bam = "{odir}/{aext}/hisat2_cuff/sample1/aligned_sorted.bam"
     log:
          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="hisat2_cuff.log")
     resources:
          memory=10
     threads:
          config["hisat2"]["threads"]
     message: """--- running hisat2"""
     shell:"""
            echo "Sample : ${input.sample} , GTF : ${input.splice_bed} , Output : ${output.outdir}" 
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/hisat2 -p {threads} \
                   -S {output.outdir}/sample1/accepted_hits.sam \
                   --dta-cufflinks \
                   --known-splicesite-infile {input.splice_bed} \  
                   -x {input.ref}/genome -f \
                   -1 {input.sample[0]} -2 {input.sample[1]} &> {output.out_log}
            cd {output.outdir}/sample1
            samtools view -bS -o aligned.bam accepted_hits.sam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam > aligned_nsort.bed
            cd -
            """

rule star:
     input:
          exe_dir=conda_dir,
          #r1,r2=inputs,
          sample=inputs,
          ref=findex("star"),
          gtf=expand("{p_gtf}/gene_exons_{{aext}}.gtf",p_gtf=gtf_p)
     output:
          outdir = "{odir}/{aext}/star",
          out_log = "{odir}/{aext}/star/sample1/run.log",
          out_bam = "{odir}/{aext}/star/sample1/aligned_sorted.bam"
     #log:
     #     logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="star.log")
     resources:
          memory=10
     threads:
          config["star"]["threads"]
     message: """--- running STAR"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/STAR --runThreadN {threads} \
                   --outFileNamePrefix {output.outdir}/sample1/sample1 \
                   --genomeDir {input.ref} \
                   --sjdbGTFfile {input.gtf} \
                   --readFilesIn {input.sample} \
                   --outFilterScoreMinOverLread 0.49                     \
                   --outFilterMatchNminOverLread 0.49 --outSAMunmapped Within \
                   --outSAMprimaryFlag OneBestScore --outSAMstrandField intronMotif \
                   --outSAMtype BAM Unsorted --quantMode GeneCounts &> {output.out_log}
            cd {output.outdir}/sample1
            mv sample1Aligned.out.bam aligned.bam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam > aligned_nsort.bed
            cd -
           """

rule hisat2_str:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=findex("hisat2"),
          splice_bed=expand("{p_gtf}/splice_sites_{{aext}}.bed",p_gtf=gtf_p)          
     output:
          outdir = "{odir}/{aext}/hisat2_str",
          out_log = "{odir}/{aext}/hisat2_str/sample1/run.log",
          out_bam = "{odir}/{aext}/hisat2_str/sample1/aligned_sorted.bam"
     log:
          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="hisat2_str.log")
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
                   -x {input.ref}/genome -f \
                   -1 {input.sample[0]} -2 {input.sample[1]} &> {output.out_log}
            cd {output.outdir}/sample1
            samtools view -bS -o aligned.bam accepted_hits.sam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam > aligned_nsort.bed
            cd -
           """  

rule star_2pass:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=findex("star"),
          gtf=expand("{p_gtf}/gene_exons_{{aext}}.gtf",p_gtf=gtf_p)
     output:
          outdir = "{odir}/{aext}/star_2pass",
          out_log = "{odir}/{aext}/star_2pass/sample1/run.log",
          out_bam = "{odir}/{aext}/star_2pass/sample1/aligned_sorted.bam"
     log:
          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="star_2pass.log")
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
            mv sample1Aligned.out.bam aligned.bam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam > aligned_nsort.bed
            cd -
           """

rule sailfish:
     input:
          exe_dir=exe_dir,
          sample=inputs,
          ref=expand("{tdir}/{{aext}}/sailfish",tdir=tfdir)
     output:
          outdir = "{odir}/{aext}/sailfish",
          out_log = "{odir}/{aext}/sailfish/sample1/run.log"
     log:
          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="sailfish.log")
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
                   -1 {input.sample[0]} -2 {input.sample[1]} \
                   -o {output.outdir}/sample1 &> {output.out_log}
           """

rule salmon:
     input:
          exe_dir=conda_dir,
          #r1,r2=inputs,
          sample=inputs,
          #ref=expand("{tdir}/{aext}/salmon",tdir=tfdir,aext=config['A_EXTN'])
          ref=expand("{tdir}/{{aext}}/salmon",tdir=tfdir)
     output:
          outdir = "{odir}/{aext}/salmon",
          out_log = "{odir}/{aext}/salmon/sample1/run.log"
     #log:
     #     logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="salmon.log")
     resources:
          memory=10
     threads:
          config["salmon"]["threads"]
     message: """--- running Salmon"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/salmon quant --numBootstraps 100 \
                   -p {threads} \
                   -i {input.ref} \
                   -l "IU" \
                   -1 {input.sample[0]} -2 {input.sample[1]} \
                   -o {output.outdir}/sample1 &> {output.out_log}
           """

rule kallisto:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          ref=expand("{tdir}/{{aext}}/kallisto/transcripts.idx",tdir=tfdir)
     output:
          outdir = "{odir}/{aext}/kallisto",
          out_log = "{odir}/{aext}/kallisto/sample1/run.log"
     log:
          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="kallisto.log")
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
                   {input.sample[0]} {input.sample[1]} &> {output.out_log}
           """
#
##soapsplice -d <2BWT findex prefix> -1 <reads_a> -2 <reads_b> -r <length of reads> -I <insert size> -o <prefix of output files> [Other Options]
#
#rule soapsplice:
#     input:
#          exe_dir=conda_dir,
#          sample=inputs,
#          ref=findex("soapsplice")
#     output:
#          outdir = "{odir}/{aext}/soapsplice",
#          out_log = "{odir}/{aext}/soapsplice/sample1/run.log"
#     log:
#          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="soapsplice.log")
#     resources:
#          memory=10
#     threads:
#          config["soapsplice"]["threads"]
#     message: """--- running soapsplice"""
#     shell:"""
#            mkdir -p {output.outdir}/sample1
#            {input.exe_dir}/soapsplice \
#                   -d {input.ref} \
#                   -1 {input.sample} -2 {input.sample} \
#                   -I 50 \
#                   -o {output.outdir}/sample1 &> {output.out_log}
#           """
#
##OUTPUT_DIR=/home/small_results/ 
##READ_FILE_END1=/home/small_1.fastq
##READ_FILE_END2=/home/small_2.fastq
##MAPSPLICE_DIR=/home/MapSplice-v2.2.0/
##REF_GENOME=/home/hg19_sequence/
##BOWTIE_INDEX=/home/hg19_findex/humanchridx_M
#
#
##python $MAPSPLICE_DIR/mapsplice.py \
##       -1 $READ_FILE_END1 \
##       -2 $READ_FILE_END2 \
##       -c $REF_GENOME \
##       -x $BOWTIE_INDEX \
##       -p 8 \
##       -o $OUTPUT_DIR 2>log.txt
#
#rule mapsplice2:
#     input:
#          exe_dir=conda_dir,
#          sample=inputs,
#          ref=get_contigs,
#          idx=findex("bowtie2")
#     output:
#          outdir = "{odir}/{aext}/mapsplice2",
#          out_log = "{odir}/{aext}/mapsplice2/sample1/run.log"
#     log:
#          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="mapsplice2.log")
#     resources:
#          memory=10
#     threads:
#          config["mapsplice"]["threads"]
#     message: """--- running mapsplice2"""
#     shell:"""
#            mkdir -p {output.outdir}/sample1
#            {input.exe_dir}/mapsplice.py -p {threads} \
#                   -1 {input.sample} -2 {input.sample} \
#                   -c {input.ref} \
#                   -x {input.idx}\
#                   -o {output.outdir}/sample1 &> {output.out_log}
#           """
#
#rule subread:
#     input:
#          exe_dir=conda_dir,
#          sample=inputs,
#          idx=findex("subread"),
#          gtf=expand("{p_gtf}/gene_exons_{{aext}}.gtf",p_gtf=gtf_p)
#     output:
#          outdir = "{odir}/{aext}/subread",
#          out_log = "{odir}/{aext}/subread/sample1/run.log"
#     log:
#          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="subread.log")
#     resources:
#          memory=10
#     threads:
#          config["subread"]["threads"]
#     message: """--- running subread"""
#     shell:"""
#            mkdir -p {output.outdir}/sample1
#            {input.exe_dir}/subread-align -T {threads} -t 0 -d 50 -D 600 \
#                   -i {input.idx}\
#                   -r {input.sample[0]} -R {input.sample[1]} \
#                   -a {input.gtf} \
#                   -o {output.outdir}/sample1/aligned.bam &> {output.out_log}
#            cd {output.outdir}/sample1
#            samtools view -h -F 4 -b -o align.bam aligned.bam
#            samtools sort -n -o aligned_nsort.bam align.bam
#            samtools sort -o aligned_sorted.bam aligned.bam
#            bedtools bamtobed -i aligned_nsort.bam > aligned_nsort.bed
#            cd -
#           """

##Paired : gsnap -d <genome> <fastq_file_1> <fastq_file_2> [<fastq_file_3> <fastq_file_4>...]
##Single  :gsnap -d <genome> --force-single-end <fastq_file_1> [<fastq_file_2>...]
##Use only $HOME/bin/gmap and $HOME/bin/gsnap; Do this before launching gsnap export LD_LIBRARY_PATH=
##cat ../annotations/gene_exons_0pcent.gtf | gtf_splicesites > gene_exons_0pcent.splicesites
##cat gene_exons_0pcent.splicesites | iit_store -o splicesites.map
##-D dir to the index path
##-d the final genome dire ( ie: genome )
##-A output to sam
## gsnap -D ../genome/0/gmap -d genome -s splicesites.map ../reads/5M/sample01_1.fasta ../reads/5M/sample01_2.fasta -A sam > gsnap/sample1.sam
##
rule gsnap:
     input:
          exe_dir = exe_dir,
          sample = inputs,
          gtf=expand("{p_gtf}/gene_exons_{{aext}}.gtf",p_gtf=gtf_p),
          idx = findex("genome")
     output:
          outdir = "{odir}/{aext}/gsnap",
          outmap = "{odir}/{aext}/gsnap/{aext}.map.iit",
          out_log = "{odir}/{aext}/gsnap/sample1/run.log"
     log:
          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="gsnap.log")
     resources:
          memory=10
     threads:
          config["gmap"]["threads"]
     message: """--- running gsnap"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            cat {input.gtf} | {input.exe_dir}/gtf_splicesites | cat - | {input.exe_dir}/iit_store -o {output.outmap}
            {input.exe_dir}/gsnap -t {threads} \
                   -D {input.idx} \
                   -d genome \
                   -s {output.outmap} \
                   {input.sample[0]} {input.sample[1]} \
                   -A sam > {output.outdir}/sample1/sample1.sam &> {output.out_log}
            cd {output.outdir}/sample1
            samtools view -bS -o aligned.bam sample1.sam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam > aligned_nsort.bed
            cd -
           """

##$ /PATH_TO_NOVOALIGN/novoalign -d /PATH_TO_PATHOGEN_GENOMES/genomes.nix -f /PATH_TO_SIMULATED_DATA/reads/simulated_set.fq -o SAM -r All > pathogen.sam
#
rule novoalign:
     input:
          exe_dir=conda_dir,
          sample=inputs,
          idx=findex("genome.nix")
     output:
          outdir = "{odir}/{aext}/novoalign",
          out_log = "{odir}/{aext}/novoalign/sample1/run.log"
     log:
          logs_path=expand("{{odir}}/{{aext}}/{tlog}",tlog="novoalign.log")
     resources:
          memory=10
     threads:
          config["novoalign"]["threads"]
     message: """--- running novoalign"""
     shell:"""
            mkdir -p {output.outdir}/sample1
            {input.exe_dir}/novoalign -p {threads} \
                   -d {input.idx} \
                   -f {input.sample[0]} {input.sample[1]} \
                   -o SAM -r All > {output.outdir}/sample1/sample1.sam &> {output.out_log}
            cd {output.outdir}/sample1
            samtools view -bS -o aligned.bam sample1.sam
            samtools view -h -F 4 -b -o align.bam aligned.bam
            samtools sort -n -o aligned_nsort.bam align.bam
            samtools sort -o aligned_sorted.bam aligned.bam
            bedtools bamtobed -i aligned_nsort.bam > aligned_nsort.bed 
            cd -
          """
