import os
import pandas as pd

# 1. Setup and Configuration
PRJNAME = config.get("PRJNAME") 
if not PRJNAME:
    raise ValueError("ERROR: Specify PRJNAME via --config")

# Load project-specific configuration from the project directory
configfile: f"reads/{PRJNAME}/config.yaml"

# 2. Include Modular Rule Files
include: "toolbox/getdata.smk" 
include: "toolbox/qc.smk" 
include: "toolbox/aligner.smk" 
include: "toolbox/vcf.smk"  
# include: "toolbox/rnaseq.smk" 

# 3. Get Functions
def get_samples(wildcards):
    """
    Trigger checkpoints base on dsource from config.yaml (default: SRA)
    local_runinfo.csv will be generated based on fastq file names under read/[project_name]
    """
    dsource = config.get("DSOURCE", "SRA").upper()
    
    # Map data source to the corresponding checkpoint
    if dsource == "SRA":
        cp_output = checkpoints.getdata_sra.get(**wildcards).output.csv 
    else:
        cp_output = checkpoints.getdata_local.get(**wildcards).output.csv 
        
    df = pd.read_csv(cp_output) 
    return df['Run'].tolist() #make list using "Run" (SRR number) column


def get_qc(wildcards):
    """
    Generates QC report targets and forces trimming.
    """
    samples = get_samples(wildcards)
    
    # Standard FastQC reports
    fastqc = expand("reads/{prj}/qc/{s}_{pair}_fastqc.html", 
                    prj=PRJNAME, s=samples, pair=[1, 2])
    
    # Trimmed fastq files (Note: Trim Galore naming uses _1_val_1 and _2_val_2)
    trim_r1 = expand("reads/{prj}/qc_trimmed/{s}_1_val_1.fq.gz", prj=PRJNAME, s=samples)
    trim_r2 = expand("reads/{prj}/qc_trimmed/{s}_2_val_2.fq.gz", prj=PRJNAME, s=samples)
    
    # Return a combined list of all targets
    return fastqc + trim_r1 + trim_r2


def get_bam(wildcards):
    """
    Generates Bam file for aligner workflow.
    """
    samples = get_samples(wildcards)
    aligner = config.get("ALIGNER", "bwa")
    return expand("reads/{prj}/bam/filtered/{s}.{aln}.filtered.bam",
                  prj=PRJNAME, s=samples, aln=aligner)


def get_vcf(wildcards):
    """
    Generates vcf.gz, merged.vcf.gz, 
    """
    samples = get_samples(wildcards)
    aligner = config.get("ALIGNER", "bwa") 
    annvcf =  f"reads/{PRJNAME}/vcf/all_samples.{aligner}.ann.vcf"
    annstats = f"reads/{PRJNAME}/vcf/all_samples.{aligner}.snpeff_stats.html"
    consensus = expand("reads/{prj}/vcf/consensus/{s}.{aln}.consensus.fa", prj=PRJNAME, s=samples, aln=aligner)
    return [annvcf, annstats] + consensus
    

def get_rnaseq(wildcards):
    """
    Generates file list for the RNA-seq workflow.
    """
    samples = get_samples(wildcards)
    return expand("reads/{prj}/counts/{s}.counts.txt", 
                  prj=PRJNAME, s=samples)



# 4. Terminal Rules
rule note:
    input:
        # Leaving this empty ensures it doesn't trigger other rules
    run:
        print("\n" + "="*50)
        print("THIS SNAKEFILE SERVICE AS TERMINAL FOR BIOINFOMATIC TOOLBOX")
        print("Please specify a target rule:")
        print("  snakemake qc      - Run FastQC")
        print("  snakemake bam     - Run Alignments")
        print("  snakemake vcf     - Run Variant Calling")
        print("  snakemake rnaseq  - Run RNA-seq analysis")
        print("="*50 + "\n")

rule qc:
    """
    QC for fastq files.
    """
    input: get_qc


rule bam:
    """
    Align fastq files with index by chosen aligner
    """
    input: get_bam

rule vcf:
    """
    Entry point for the Variant Calling workflow.
    """
    input: get_vcf

rule rnaseq:
    """
    Entry point for the RNA-seq workflow.
    """
    input: get_rnaseq








