from functools import reduce
import pandas as pd

configfile: 'config.yaml'

################################################################################
# Globals                                                                      #
################################################################################
exe_dir = config['EXE_DIR']
conda_dir = config['CONDA_DIR']
base_dir = config['BASE']
org = config['ORG']
samples = pd.read_csv(config['sample_file'], sep=' ')
is_paired = "fasta2" in samples.columns
etype = [ "alt_genome_annotations" , "alt_coverage" , "alt_reference"]
exp_type = etype[0] 

################################################################################
# Tools used in the
################################################################################

aligners = ["tophat2", "hisat2", "star" , "gsnap" , "novoalign" , "soapsplice" , "mapsplice", "kallisto", "sailfish", "salmon" ]
assemblers = ["cufflinks" , "stringtie" , "RSEM" , "kallisto", "sailfish" , "salmon" , "htseq" ,"express" ]

################################################################################
# Set Input/ Output directory paths
################################################################################

reads_dir = base_dir + "/" + org + "/reads/" + config['TRUE_COV']
#index_base = base_dir + "/" + org + "/genome/0/"
#index_base = base_dir + "/" + org + "/genome/0/"
results_base = base_dir + "/" + org + "/results/" + exp_type
logs_base = base_dir  + "/" + org + "/logs/" + exp_type
ann_results_dir = [ results_base +"/"+x for x in config['A_EXTN'] ]
ann_logs_dir = [ logs_base +"/"+x for x in config['A_EXTN'] ]

# Remove variables #
#gt_logs_dir = [ logs_dir[2]+"/"+x for x in config['GT_EXTN'] ]
#cov_logs_dir = [ logs_dir[1]+"/"+x for x in config['COV_EXTN'] ]

## Make one directory variable listing all directories in the experiment
dirs= [j for i in [ann_results_dir,ann_logs_dir] for j in i]

################################################################################
# Functions                                                                    #
################################################################################

def get_samples():
    """Returns list of all samples."""
    return list(samples["sample"].unique())


def get_samples_with_rep():
    """Returns list of all combined condition/sample identifiers."""
    return list((samples["sample"] + "." + samples["condition"]).unique())


def get_samples_by_condition(condition):
    """Returns lanes for given sample."""
    subset = samples.loc[samples["condition"] == condition]
    return list(subset["sample"].unique())

################################################################################
# Rules                                                                        #
################################################################################

##### setup report #####

report: "report/workflow.rst"

########################

rule all:
    input:
        dirs,
        expand("{rbase}/{aext}/{aligners}/run.log",rbase=results_base,aext=config['A_EXTN'],aligners=aligners)
        #"counts/merged.log2.txt",
        #"qc/multiqc_report.html"

rule mkdirs:
    output:
         dirs
    shell:
         "mkdir -p "+' '.join(dirs)

include: "rules/alignment.smk"
#include: "rules/assemble.smk"
#include: "rules/alignment_qc.smk"
#include: "rules/counts.smk"
