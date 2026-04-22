localrules: fastqc, trim_reads
rule fastqc:
    input:
        r1 = f"reads/{PRJNAME}/{{sample}}_1.fastq",
        r2 = f"reads/{PRJNAME}/{{sample}}_2.fastq"
    output:
        html1 = f"reads/{PRJNAME}/qc/{{sample}}_1_fastqc.html",
        html2 = f"reads/{PRJNAME}/qc/{{sample}}_2_fastqc.html",
        zip1  = f"reads/{PRJNAME}/qc/{{sample}}_1_fastqc.zip",
        zip2  = f"reads/{PRJNAME}/qc/{{sample}}_2_fastqc.zip"
    log: f"reads/{PRJNAME}/logs/qc/fastqc/{{sample}}.log"
    conda: "../env/qc.yaml"
    shell:
        """
        # -s checks if the file exists and has size > 0
        if [ -s "{input.r2}" ]; then
            fastqc {input.r1} {input.r2} --outdir reads/{PRJNAME}/qc/ > {log} 2>&1
        else
            # Single-end: only run R1
            fastqc {input.r1} --outdir reads/{PRJNAME}/qc/ > {log} 2>&1
            
            # Create dummy outputs for R2 to satisfy Snakemake
            touch {output.html2} {output.zip2}
            echo "R2 is empty dummy, FastQC skipped R2." >> {log}
        fi
        """

rule multiqc:
    input:
        # Forces the rule to wait until all individual FastQC reports exist
        lambda wildcards: expand("reads/{prj}/qc/{s}_{pair}_fastqc.html", 
                                 prj=PRJNAME, s=get_samples(wildcards), pair=[1, 2])
    output:
        f"reads/{PRJNAME}/qc/multiqc_report.html"
    params:
        qc_dir = f"reads/{PRJNAME}/qc/"
    conda: "../env/multiqc.yaml"
    shell:
        # --force overwrites old reports, so the program wont stop if old one exist
        "multiqc {params.qc_dir} -o {params.qc_dir} -n multiqc_report.html --force"

rule trim_reads:
    """
    Handles quality trimming (tail cutting) at Phred 30.
    """
    input:
        r1 = f"reads/{PRJNAME}/{{sample}}_1.fastq",
        r2 = f"reads/{PRJNAME}/{{sample}}_2.fastq"
    output:
        # Trim Galore naming convention for paired: _val_1 and _val_2
        r1_p = f"reads/{PRJNAME}/qc_trimmed/{{sample}}_1_val_1.fq.gz",
        r2_p = f"reads/{PRJNAME}/qc_trimmed/{{sample}}_2_val_2.fq.gz"
    log: f"reads/{PRJNAME}/logs/qc/trim_reads/{{sample}}.log"
    conda: "../env/qc.yaml"
    shell:
        """
        if [ -s "{input.r2}" ]; then
            # Real paired-end data
            trim_galore --paired --gzip -q 30 --output_dir reads/{PRJNAME}/qc_trimmed/ {input.r1} {input.r2} > {log} 2>&1
        else
            # Single-end data: Trim Galore creates a different filename for SE
            trim_galore --gzip -q 30 --output_dir reads/{PRJNAME}/qc_trimmed/ {input.r1} > {log} 2>&1
            
            # 1. Rename the SE output to match expected file name
            mv reads/{PRJNAME}/qc_trimmed/{wildcards.sample}_1_trimmed.fq.gz {output.r1_p}
            
            # 2. Create a dummy gzipped file for R2 so downstream rules don't crash
            echo "" | gzip > {output.r2_p}
            echo "Single-end mode: Renamed R1 and created dummy R2.gz" >> {log}
        fi
        """