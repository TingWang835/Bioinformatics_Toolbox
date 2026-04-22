# 1. Checkpoints
localrules: getdata_sra, getdata_local, download_refs, samtools_faidx, download_fastq, prepare_reference

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
    output: f"reads/{PRJNAME}/local_runinfo.csv"
    shell:
        """
        echo "Run,LibraryLayout" > {output.csv}
        for f in reads/{PRJNAME}/*_1.fast*; do
            base=$(basename $f | sed 's/_1.fast.*//')

            # Check if _2.fast* exists, determine library layout (SE/PE)
            if [ -s "reads/{PRJNAME}/${{base}}_2.fast*" ]; then
                echo "$base,PAIRED" >> {output.csv}
            else
                echo "$base,SINGLE" >> {output.csv}
            fi
        
        done
        """        


# 2. rules
rule download_refs:
    output:
        fasta = f"refs/{config['REFNAME']}/{config['ACC']}.fa",
        gff = f"refs/{config['REFNAME']}/{config['ACC']}.gff"
    params: acc = f"{config['ACC']}"
    conda: "../env/getdata.yaml"
    log: f"reads/{PRJNAME}/logs/getdata/download_refs_{config['ACC']}.log"
    shell:
        """
        search_result=$(esearch -db nuccore -query {params.acc})
        
        # 1. Fetch FASTA (for aligner and vcf annotation)
        echo "$search_result" | efetch -format fasta > {output.fasta}
        
        # 2. Fetch GFF3 (for vcf annotation)
        echo "$search_result" | efetch -format gff3 > {output.gff}
        """

## this is a unfinished function used in adding metadata to gff files, which can later enhance snpeff_build to create better annotation db
# rule augment_gff:
#     input: f"refs/{config['REFNAME']}/{config['ACC']}.gff"
#     output: f"refs/{config['REFNAME']}/{config['ACC']}.enhanced.gff"
#     shell:
#         # Example: Adding a project tag to every gene line
#         "sed 's/gene_id=/project={PRJNAME};gene_id=/g' {input} > {output}"

rule download_fastq:
    output:
        r1 = f"reads/{PRJNAME}/{{sample}}_1.fastq",
        r2 = f"reads/{PRJNAME}/{{sample}}_2.fastq"
    conda: "../env/getdata.yaml"
    log: f"reads/{PRJNAME}/logs/getdata/download_fastq/{{sample}}.log"
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


rule samtools_faidx:
    input:
        ref = f"refs/{config['REFNAME']}/{config['ACC']}.fa"
    output:
        fai = f"refs/{config['REFNAME']}/{config['ACC']}.fa.fai"
    conda: "../env/aligner.yaml"
    log: f"reads/{PRJNAME}/logs/getdata/samtools_faidx_{{config['ACC']}}.log"
    shell:
        "samtools faidx {input.ref} > {log} 2>&1"

rule prepare_reference: # shortcut to run download_refs and samtools_faidx
    input:
        f"refs/{config['REFNAME']}/{config['ACC']}.fa",
        f"refs/{config['REFNAME']}/{config['ACC']}.fa.fai"