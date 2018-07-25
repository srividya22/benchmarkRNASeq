import os

# Benchmarking Assemblers 


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

bam_tools = [ "tophat2" ,"hisat2_cuff", "hisat2_str", "star" , "star_2pass" ]
#print(get_gtf)
#print(get_outdir)
#print(get_logsdir)

rule cufflinks:
     input: 
          exe_dir=exe_dir,
          sample=expand("{rbase}/{aext}/aligner/{tool}/aligned_sorted.bam",rbase=results_base,aext=config['A_EXTN'],tool=bam_tools)
          gtf=expand("{p_gtf}/gene_exons_{aext}.gtf",p_gtf=gtf_p,aext=config['A_EXTN'])          
     output:
          outdir = expand("results_base/{aext}/assembly/cufflinks",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/assembly/cufflinks/run.log",odir=results_base,aext=config['A_EXTN']),
          out_gtf = expand("{odir}/{aext}/assembly/cufflinks/sample1/transcripts.gtf",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/assembly/cufflinks.log",odir=results_base,aext=config['A_EXTN'])
     resources:
          memory=10
     threads:
          config["cufflinks"]["threads"]
     message: """--- running Cufflinkss"""
     shell:
            """ 
            echo "Sample : ${{input.sample}} , GTF : ${{input.gtf}} , Output : ${{output.outdir}}" 
            mkdir -p {output.outdir}
            exe_dir/cufflinks \
                   -o {output.outdir}/sample1 \
                   -G {input.gtf} \
                   -p {threads} \
                      {input.sample} &> {output.out_log}
            """

rule stringtie:
     input:
          exe_dir=exe_dir,
          sample=expand("{rbase}/{aext}/aligner/{tool}/aligned_sorted.bam",rbase=results_base,aext=config['A_EXTN'],tool=bam_tools)
          gtf=expand("{p_gtf}/gene_exons_{aext}.gtf",p_gtf=gtf_p,aext=config['A_EXTN'])
     output:
          outdir = expand("results_base/{aext}/assembly/stringtie",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/assembly/stringtie/run.log",odir=results_base,aext=config['A_EXTN']),
          out_gtf = expand("{odir}/{aext}/assembly/stringtie/sample1/transcripts.gtf",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/assembly/stringtie.log",odir=results_base,aext=config['A_EXTN'])
     resources:
          memory=10
     threads:
          config["stringtie"]["threads"]
     message: """--- running Stringtie"""
     shell:
            """
            echo "Sample : ${{input.sample}} , GTF : ${{input.gtf}} , Output : ${{output.outdir}}"
            mkdir -p {output.outdir}/sample1
            conda_dir/stringtie \
                   -o {output.outdir}/sample1/transcripts.gtf \
                   -G {input.gtf} \
                   -p {threads} \
                   -A {output.outdir}/sample1/genes_fpkm.tracking" \
                   -b {output.outdir}/sample1 \
                      {input.sample} &> {output.out_log}
            """

rule featurecounts:
     input:
          exe_dir=exe_dir,
          sample=expand("{rbase}/{aext}/aligner/{tool}/aligned_sorted.bam",rbase=results_base,aext=config['A_EXTN'],tool=bam_tools)
          gtf=expand("{p_gtf}/gene_exons_{aext}.gtf",p_gtf=gtf_p,aext=config['A_EXTN'])
     output:
          outdir = expand("results_base/{aext}/assembly/featurecounts",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/assembly/featurecounts/run.log",odir=results_base,aext=config['A_EXTN']),
          out_gtf = expand("{odir}/{aext}/assembly/featurecounts/sample1/transcripts_counts.txt",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/assembly/featurecounts.log",odir=results_base,aext=config['A_EXTN'])
     resources:
          memory=10
     threads:
          config["featurecounts"]["threads"]
     message: """--- running featurecounts"""
     shell:
            """
            echo "Sample : ${{input.sample}} , GTF : ${{input.gtf}} , Output : ${{output.outdir}}"
            mkdir -p {output.outdir}
            conda_dir/featureCounts \
                   -T {threads} \
                   -Q 10 \
                   -p -M --fraction -J -g transcript_id \
                   -a {input.gtf} \
                   -o {output.outdir}/sample1/transcripts_counts.txt \
                      {input.sample} &> {output.out_log}
            """

rule RSEM:
     input:
          exe_dir=exe_dir,
          sample=expand("{rbase}/{aext}/aligner/{tool}/aligned_sorted.bam",rbase=results_base,aext=config['A_EXTN'],tool=bam_tools)
          gtf=expand("{p_gtf}/gene_exons_{aext}.gtf",p_gtf=gtf_p,aext=config['A_EXTN'])
     output:
          outdir = expand("results_base/{aext}/assembly/rsem",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/assembly/rsem/run.log",odir=results_base,aext=config['A_EXTN']),
          out_gtf = expand("{odir}/{aext}/assembly/rsem/transcripts.gtf",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/assembly/rsem.log",odir=results_base,aext=config['A_EXTN'])
     resources:
          memory=10
     threads:
          config["rsem"]["threads"]
     message: """--- running Cufflinkss"""
     shell:
            """
            echo "Sample : ${{input.sample}} , GTF : ${{input.gtf}} , Output : ${{input.outdir}}"
            mkdir -p {output.outdir}
            conda_dir/RSEM \
                   -o {output.outdir}/sample1 \
                   -G {input.gtf} \
                   -p {threads} \
                      {input.sample} &> {output.out_log}
            """

rule express:
     input:
          exe_dir=exe_dir,
          sample=expand("{rbase}/{aext}/aligner/{tool}/aligned_sorted.bam",rbase=results_base,aext=config['A_EXTN'],tool=bam_tools)
          gtf=expand("{p_gtf}/gene_exons_{aext}.gtf",p_gtf=gtf_p,aext=config['A_EXTN'])
     output:
          outdir = expand("results_base/{aext}/assembly/express",odir=results_base,aext=config['A_EXTN']),
          out_log = expand("{odir}/{aext}/assembly/express/run.log",odir=results_base,aext=config['A_EXTN']),
          out_gtf = expand("{odir}/{aext}/assembly/express/transcripts.gtf",odir=results_base,aext=config['A_EXTN'])
     log:
          logs_path=expand("{odir}/{aext}/assembly/express.log",odir=results_base,aext=config['A_EXTN'])
     resources:
          memory=10
     threads:
          config["express"]["threads"]
     message: """--- running Cufflinkss"""
     shell:
            """
            echo "Sample : ${{input.sample}} , GTF : ${{input.gtf}} , Output : ${{input.outdir}}"
            mkdir -p {output.outdir}
            conda_dir/express \
                   -o {output.outdir}/sample1 \
                   -G {input.gtf} \
                   -p {threads} \
                      {input.sample} &> {output.out_log}
