localrules: getdata_sra, getdata_local, download_ncbi_ref, samtools_faidx, download_fastq

# =============================================================================
# Check Points
# =============================================================================

checkpoint getdata_sra:
    conda: "../env/getdata.yaml"
    output: csv = f"{READS_DIR}/sra_runinfo.csv"
    shell:
        """
        if [ -f "{output.csv}" ]; then
            echo "sra_runinfo already exists, Loading {output.csv}"  
        else
            esearch -db sra -query {config[PRJNUMBER]} | efetch -format runinfo > {output.csv}
        fi
        """
        
checkpoint getdata_local:
    output: csv = f"{READS_DIR}/local_runinfo.csv"
    shell:
        """
        echo "Run,LibraryLayout" > {output.csv}
        for f in {READS_DIR}/*_1.fast*; do
            base=$(basename $f | sed 's/_1.fast.*//')
            if [ -s "{READS_DIR}/${{base}}_2.fast*" ]; then
                echo "$base,PAIRED" >> {output.csv}
            else
                echo "$base,SINGLE" >> {output.csv}
            fi
        done
        """

# =============================================================================
# Helper Functions
# =============================================================================
def get_ref_source(wildcards):
    """Selects the correct FASTA, GTF, and GFF file paths based on the configured provider."""
    if ref_src == "ensembl":
        return {
            "ref": f"{REFS_DIR}/{species_cap}.{assembly}.dna.toplevel.fa",
            "gtf": f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gtf",
            "gff": f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gff3"
        }
    else:
        return {
            "ref": f"{REFS_DIR}/{acc}.fa",
            "gtf": f"{REFS_DIR}/{acc}.gtf",
            "gff": f"{REFS_DIR}/{acc}.gff"
        }

# =============================================================================
# Rules
# =============================================================================
rule download_ncbi_ref:
    output:
        fasta = f"{REFS_DIR}/{acc}.fa", 
        gff   = f"{REFS_DIR}/{acc}.gff" 
    params:
        search_cmd = (
            f"esearch -db assembly -query {acc}| elink -target nuccore"
            if acc.startswith(("GCF_", "GCA_")) else
            f"esearch -db nuccore -query {acc}"
        )
    conda: "../env/getdata.yaml" 
    log: f"{LOG_DIR}/getdata/download_ncbi_ref_{acc}.log" 
    shell:
        """
        search_result=$({params.search_cmd})
        echo "$search_result" | efetch -format fasta > {output.fasta} 2> {log}
        echo "$search_result" | efetch -format gff3 > {output.gff} 2>> {log}
        """

rule download_ensembl_refs:
    output:
        fa   = f"{REFS_DIR}/{species_cap}.{assembly}.dna.toplevel.fa",
        gff3 = f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gff3"
    log: f"{LOG_DIR}/getdata/download_ensembl_refs.log"
    shell:
        """
        wget -O {output.fa}.gz ftp://ftp.ensembl.org/pub/release-{release}/fasta/{species}/dna/{species_cap}.{assembly}.dna.toplevel.fa.gz > {log} 2>&1
        wget -O {output.gff3}.gz ftp://ftp.ensembl.org/pub/release-{release}/gff3/{species}/{species_cap}.{assembly}.{release}.gff3.gz >> {log} 2>&1
        gunzip -f {output.fa}.gz
        gunzip -f {output.gff3}.gz
        """

rule download_fastq:
    output:
        r1 = f"{READS_DIR}/{{sample}}_1.fastq",
        r2 = f"{READS_DIR}/{{sample}}_2.fastq"
    conda: "../env/getdata.yaml"
    log: f"{LOG_DIR}/getdata/download_fastq/{{sample}}.log"
    params:
        dsource = config.get("DATASOURCE", "SRA").upper(),
        n = config.get("N", 10000)
    shell:
        """
        if [ "{params.dsource}" == "SRA" ]; then
            fastq-dump -X {params.n} --split-files --outdir {READS_DIR} {wildcards.sample}
            if [ ! -f "{output.r2}" ]; then touch {output.r2}; fi
        elif [ "{params.dsource}" == "LOCAL" ]; then
            if [ ! -f "{output.r1}" ]; then 
                echo "Error: Local R1 file {output.r1} missing!";
                exit 1
            fi
            if [ ! -f "{output.r2}" ]; then touch {output.r2}; fi
        else
            echo "Error: Invalid DATASOURCE '{params.dsource}'. Use SRA or LOCAL."
            exit 1
        fi
        """

rule gff_to_gtf:
    input:
        gff3 = (
            f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gff3"
            if ref_src == "ensembl" else f"{REFS_DIR}/{acc}.gff"
        )
    output:
        gtf = (
            f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gtf"
            if ref_src == "ensembl" else f"{REFS_DIR}/{acc}.gtf"
        )
    log: f"{LOG_DIR}/getdata/gff_to_gtf.log"
    conda: "../env/getdata.yaml"
    shell:
        """
        gffread {input.gff3} -T -o {output.gtf} > {log} 2>&1
        """

rule samtools_faidx:
    input:
        unpack(get_ref_source)
    output:
        fai = (
            f"{REFS_DIR}/{species_cap}.{assembly}.dna.toplevel.fa.fai"
            if ref_src == "ensembl" else f"{REFS_DIR}/{acc}.fa.fai"
        )
    log: f"{LOG_DIR}/getdata/samtools_faidx.log"
    conda: "../env/dna_aligner.yaml"
    shell:
        """
        samtools faidx {input.ref} > {log} 2>&1
        """