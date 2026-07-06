localrules: getdata_sra, getdata_local, download_ncbi_ref, samtools_faidx, download_fastq

# =============================================================================
# Check Points
# =============================================================================

checkpoint getdata_sra:
    """
    Fetches runinfo details via SRA project numbers.
    """
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
    """
    Generate runinfo details from local filenames.
    """
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
# Rules
# =============================================================================
rule download_ncbi_ref:
    output:
        fasta = f"{REFS_DIR}/{acc}.fa", 
        gff   = f"{REFS_DIR}/{acc}.gff3" 
    conda: "../env/getdata.yaml" 
    log: f"{LOG_DIR}/getdata/download_ncbi_ref.log" 
    shell:
        """
        # DO NOT REMOVE mkdir
        mkdir -p {REFS_DIR} 
        exec 2> {log}

        TMP_ZIP="{REFS_DIR}/{acc}_dataset.zip"

        # Download zip (include fa and gff) via NCBI datasets API CLI
        datasets download genome accession {acc} \
            --include genome,gff3 \
            --filename $TMP_ZIP

        # Unzip package
        unzip -q -o $TMP_ZIP -d {REFS_DIR}/tmp_out

        # Move fa and gff to path
        mv {REFS_DIR}/tmp_out/ncbi_dataset/data/{acc}/*_genomic.fna {output.fasta}
        mv {REFS_DIR}/tmp_out/ncbi_dataset/data/{acc}/genomic.gff {output.gff}

        # Purge temporary packages
        rm -rf $TMP_ZIP {REFS_DIR}/tmp_out
        """

rule download_ensembl_refs:
    output:
        fasta = f"{REFS_DIR}/{species_cap}.{assembly}.{release}.fa",
        gff   = f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gff3"
    log:
        f"{LOG_DIR}/getdata/download_ensembl_ref.log"
    conda: "../env/getdata.yaml"
    params:
        url_base = f"ftp://ftp.ensembl.org/pub/release-{release}"
    shell:
        """
        exec 2> {log}

        # 1. Stream, decompress, and rename the Genomic FASTA file
        wget -qO- {params.url_base}/fasta/{species_low}/dna/{species_cap}.{assembly}.dna.toplevel.fa.gz | \
            gunzip -c > {output.fasta}

        # 2. Stream, decompress, and rename the GFF3 Annotation file
        wget -qO- {params.url_base}/gff3/{species_low}/{species_cap}.{assembly}.{release}.gff3.gz | \
            gunzip -c > {output.gff}
        """


rule download_fastq:
    output:
        r1 = f"{READS_DIR}/{{sample}}_1.fastq",
        r2 = f"{READS_DIR}/{{sample}}_2.fastq"
    conda: "../env/getdata.yaml"
    log: f"{LOG_DIR}/getdata/download_fastq/{{sample}}.log"
    params:
        dsource = config.get("DATASOURCE", "SRA").upper(),
    shell:
        """
        if [ "{params.dsource}" == "SRA" ]; then
            fastq-dump --split-files --outdir {READS_DIR} {wildcards.sample}
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
    """
    Locally converts GFF3  to GTF format via gffread.
    """
    input:
        gff = f"{REFS_DIR}/{{ref_file}}.gff3" # Uses standardized REFS_DIR variable
    output:
        gtf = f"{REFS_DIR}/{{ref_file}}.gtf"
    log: f"{LOG_DIR}/getdata/gff_to_gtf_{{ref_file}}.log" # Named uniquely per target
    conda: "../env/getdata.yaml"
    shell:
        """
        gffread {input.gff} -T -o {output.gtf} > {log} 2>&1
        """

rule samtools_faidx:
    """
    Generates .fai index tracks on reference fasta structures.
    """
    input:
        fasta = f"{REFS_DIR}/{{ref_file}}.fa"
    output:
        fai = f"{REFS_DIR}/{{ref_file}}.fa.fai"
    log: f"{LOG_DIR}/getdata/samtools_faidx_{{ref_file}}.log" 
    conda: "../env/dna_aligner.yaml" 
    shell:
        """
        samtools faidx {input.fasta} > {log} 2>&1
        """


# import sys
# sys.path.append("toolbox/scripts")
# from uscs_chromosomes import CHROM_DICTS

# rule uscs_translator: 
# """this is a translator for generating a fasta file with USCS chrom naming system from ncbi or ensemble ref"""
#     input:
#         unpack(get_refs)
#     output:
#         uscs_fasta = f"{REFS_DIR}/uscs_chroms.fa"
#     params:
#         # Evaluates cleanly using your isolated helper script
#         mapping = get_ucsc_mapping(species_low)
#     log:
#         f"{LOG_DIR}/refs/uscs_translator.log"
#     run:
#         with open(input.fasta, "r") as infile, open(output.uscs_fasta, "w") as outfile:
#             for line in infile:
#                 if line.startswith(">"):
#                     # Extract the current header name (e.g., ">I" or ">I description")
#                     original_header = line.strip().split()[0].replace(">", "")
                    
#                     # Look up translation; fall back to original header if missing
#                     new_header = params.mapping.get(original_header, original_header)
                    
#                     # Reconstruct the header line
#                     outfile.write(f">{new_header}\n")
#                 else:
#                     outfile.write(line)