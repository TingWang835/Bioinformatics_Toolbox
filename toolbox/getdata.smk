# 1. Checkpoints
checkpoint getdata_sra:
    conda: "../env/getdata.yaml"
    output: csv = f"reads/{PRJNAME}/sra_runinfo.csv"
    shell:
        """
        if [ -f "{output.csv}" ]; then
            echo "sra_runinfo already exist, Loading {output.csv}"  
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
# 2. rules
rule download_refs:
    output:
        fasta=f"refs/{{refname}}/{{acc}}.fa",
        gff=f"refs/{{refname}}/{{acc}}.gff"
    log: f"refs/{{refname}}/logs/rule_dowload_refs/{{acc}}.log"
    conda: "../env/getdata.yaml"
    shell:
        """
        esearch -db nuccore -query {wildcards.acc} | efetch -format fasta > {output.fasta}
        esearch -db nuccore -query {wildcards.acc} | efetch -format gff > {output.gff}
        """

rule download_fastq:
    output:
        r1 = f"reads/{PRJNAME}/{{sample}}_1.fastq",
        r2 = f"reads/{PRJNAME}/{{sample}}_2.fastq"
    log: f"reads/{PRJNAME}/logs/rule_download_fastq/{{sample}}.log"
    conda: "../env/getdata.yaml"
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