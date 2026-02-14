# how to run: snakemake --cores 8 --use-conda --config PRJNAME=test DSOURCE=LOCAL -p
import os
import pandas as pd

# 1. Setup and Config
PRJNAME = config.get("PRJNAME")
if not PRJNAME:
    raise ValueError("ERROR: Specify PRJNAME via --config")

configfile: f"reads/{PRJNAME}/config.yaml"

# 2. Checkpoints
checkpoint getdata_sra:
    conda: "env/getdata.yaml"
    output: csv = f"reads/{PRJNAME}/sra_runinfo.csv"
    shell:
        """
        if [ -f "{output.csv}" ]; then
            echo "Loading existing {output.csv}"
        else
            esearch -db sra -query {config[PRJNUMBER]} | efetch -format runinfo > {output.csv}
        fi
        """

checkpoint getdata_local:
    output: csv = f"reads/{PRJNAME}/local_runinfo.csv"
    shell:
        """
        echo "Run,LibraryLayout" > {output.csv}
        for f in reads/{PRJNAME}/*_1.fastq; do
            base=$(basename $f _1.fastq)
            # Check if _2.fastq exists AND is not an empty dummy file
            if [ -s "reads/{PRJNAME}/${{base}}_2.fastq" ]; then
                echo "$base,PAIRED" >> {output.csv}
            else
                echo "$base,SINGLE" >> {output.csv}
            fi
        done
        """

# 3. Targets
def target_tag(wildcards):
    return f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}"

def target_full_alignment(wildcards):
    dsource = config.get("DSOURCE", "SRA").upper()
    aligner = config.get("ALIGNER", "bwa")
    
    if dsource == "SRA":
        checkpoint_output = checkpoints.getdata_sra.get(**wildcards).output[0]
    else:
        checkpoint_output = checkpoints.getdata_local.get(**wildcards).output[0]
        
    df = pd.read_csv(checkpoint_output)
    return [f"reads/{PRJNAME}/bam/{s}.{aligner}.bam" for s in df['Run']]

# 4. Rules
rule all:
    input: target_full_alignment

rule download_refs:
    output:
        fasta=f"refs/{{refname}}/{{acc}}.fa",
        gff=f"refs/{{refname}}/{{acc}}.gff"
    conda: "env/getdata.yaml"
    shell:
        """
        esearch -db nuccore -query {wildcards.acc} | efetch -format fasta > {output.fasta}
        esearch -db nuccore -query {wildcards.acc} | efetch -format gff > {output.gff}
        """

rule download_fastq:
    output:
        r1 = f"reads/{PRJNAME}/{{sample}}_1.fastq",
        r2 = f"reads/{PRJNAME}/{{sample}}_2.fastq"
    conda: "env/getdata.yaml"
    params:
        dsource = config.get("DSOURCE", "SRA").upper(),
        n = config.get("N", 10000)
    shell:
        """
        if [ "{params.dsource}" == "SRA" ]; then
            fastq-dump -X {params.n} --split-files --outdir reads/{PRJNAME} {wildcards.sample}
            # Create dummy R2 only if a real one does not exist
            if [ ! -f "{output.r2}" ]; then touch {output.r2}; fi

        elif [ "{params.dsource}" == "LOCAL" ]; then
            if [ ! -f "{output.r1}" ]; then 
                echo "Error: Local R1 file {output.r1} missing!"; exit 1
            fi
            # Create dummy R2 only if a real one does not exist
            if [ ! -f "{output.r2}" ]; then 
                touch {output.r2}
            fi

        else
            echo "Error: Invalid DSOURCE '{params.dsource}'. Use SRA or LOCAL."
            exit 1
        fi
        """
rule bwa_index:
    input:
        ref = f"refs/{config['REFNAME']}/{config['ACC']}.fa"
    output:
        # Using multiext to track all 5 files BWA creates
        multiext(f"refs/{config['REFNAME']}/bwa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    conda: "env/aligner.yaml"
    params:
        prefix = f"refs/{config['REFNAME']}/bwa"
    shell:
        "bwa index {input.ref} -p {params.prefix}"

rule bwa:
    input:
        r1 = f"reads/{PRJNAME}/{{sample}}_1.fastq",
        r2 = f"reads/{PRJNAME}/{{sample}}_2.fastq",
        ref = f"refs/{config.get('REFNAME', 'tmp')}/{config.get('ACC', 'tmp')}.fa",
        idx = multiext(f"refs/{config['REFNAME']}/bwa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam = temp(f"reads/{PRJNAME}/bam/{{sample}}.bwa.bam"),
        bai = (f"reads/{PRJNAME}/bam/{{sample}}.bwa.bam.bai"),
    conda: "env/aligner.yaml"
    params:
        tag = target_tag,
        refindex = f"refs/{config.get('REFNAME')}/bwa",
    shell:
        """
        # Determine if r2 is actually a data file or just a touched empty file
        if [ -s "{input.r2}" ]; then
            bwa mem -R '{params.tag}' {params.refindex} {input.r1} {input.r2} | samtools sort -o {output.bam}
        else
            bwa mem -R '{params.tag}' {params.refindex} {input.r1} | samtools sort -o {output.bam}
        fi
        samtools index {output.bam}
        """

