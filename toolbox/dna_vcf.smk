localrules: vcfcall, vcfmerge, vcfnorm, vcfconsensus, snpeff_build, vcfannotation, vcf_interactive_query, vcf_filter_by_query

# =============================================================================
# Rules
# =============================================================================

rule vcfcall:
    """
    Performs variant calling using freebayes and sorts/compresses the stream via bcftools.
    """
    input:
        unpack(get_refs),
        bam   = f"{READS_DIR}/bam/filtered/{{sample}}.{{aligner}}.filtered.bam"
    output:
        vcf = temp(f"{READS_DIR}/vcf/raw/{{sample}}.{{aligner}}.vcf.gz"),
        tbi = temp(f"{READS_DIR}/vcf/raw/{{sample}}.{{aligner}}.vcf.gz.tbi")
    params:
        ploidy = config.get("PLOIDY", 2)
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/vcfcall/{{sample}}.{{aligner}}.log"
    shell:
        """
        freebayes -f {input.fasta} --ploidy {params.ploidy} {input.bam} | bcftools view -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """

rule vcfmerge:
    """
    Merges individual sample VCF files into a multi-sample cohort VCF dataset.
    """
    input:
        vcfs = lambda w: expand(f"{READS_DIR}/vcf/raw/{{sample}}.{w.aligner}.vcf.gz", sample=get_runinfo(w)),
        tbis = lambda w: expand(f"{READS_DIR}/vcf/raw/{{sample}}.{w.aligner}.vcf.gz.tbi", sample=get_runinfo(w))
    output:
        merged_vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.vcf.gz"
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/vcfmerge/all_sample.{{aligner}}.merged.log"
    shell:
        "bcftools merge {input.vcfs} -Oz -o {output.merged_vcf}"

rule vcfnorm:
    """
    Left-aligns and normalizes indels against the reference sequence assembly.
    """
    input:
        unpack(get_refs),
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.vcf.gz"
    output:
        vcf_gz = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz.tbi"
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/vcfnorm/all_samples.{{aligner}}.merged.norm.log"
    shell:
        """
        bcftools norm -f {input.fasta} -m -any {input.vcf} -Oz -o {output.vcf_gz}
        bcftools index -t {output.vcf_gz}
        """

rule vcfconsensus:
    """
    Applies genomic variants onto the reference assembly to extract sample-specific consensus sequences.
    """
    input:
        unpack(get_refs),
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz.tbi"
    output:
        fasta = f"{READS_DIR}/vcf/consensus/{{sample}}.{{aligner}}.consensus.fa"
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/vcfconsensus/{{sample}}.{{aligner}}.consensus.log"
    shell:
        "bcftools consensus -f {input.fasta} -s {wildcards.sample} {input.vcf} > {output.fasta}"

# =============================================================================
# Functional Annotation & Queries
# =============================================================================

rule snpeff_build:
    """
    Builds a local database instance inside snpEff structure environments.
    """
    input:
        unpack(get_refs)
    output:
        db = f"{DATABASES_DIR}/snpEffectPredictor.bin"
    params:
        data_dir = lambda w: os.path.dirname(DATABASES_DIR.rstrip("/")),
        # Drops the path prefix to extract the exact folder name Snakemake expects
        genome = lambda w: os.path.basename(DATABASES_DIR.rstrip("/"))
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/snpeff_build.log"
    shell:
        """
        # Step 1: create the precise nested genome subdirectory
        mkdir -p {params.data_dir}/{params.genome}
        cp {input.fasta} {params.data_dir}/{params.genome}/sequences.fa
        cp {input.gff} {params.data_dir}/{params.genome}/genes.gff

        # Step 2: create local snpeff.config inside the genome directory
        ABS_PATH=$(pwd)
        echo "{params.genome}.genome : {params.genome}" > {params.data_dir}/{params.genome}/snpEff.config

        # Step 3: building DB pointing to the parent data directory
        snpEff build -gff3 -v {params.genome} \
            -c $ABS_PATH/{params.data_dir}/{params.genome}/snpEff.config \
            -dataDir $ABS_PATH/{params.data_dir} \
            -noCheckCds -noCheckProtein -nodownload > {log} 2>&1

        # Step 4: cleaning up and organizing for vcfannotation compatibility
        if [ -f {output.db} ]; then
            rm {params.data_dir}/{params.genome}/genes.gff
        else
            exit 1
        fi
        """

rule vcfannotation:
    """
    Annotates localized functional variant alerts using local snpEff predictor profiles.
    """
    input:
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz",
        db = f"{DATABASES_DIR}/snpEffectPredictor.bin"
    output:
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.ann.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{{aligner}}.ann.vcf.gz.tbi",
        stats = f"{READS_DIR}/vcf/all_samples.{{aligner}}.snpeff_stats.html"
    params:
        genome = lambda w: os.path.basename(DATABASES_DIR.rstrip("/")),
        ram = config.get("RAM", "4g"),
        db_dir = lambda w: os.path.dirname(DATABASES_DIR.rstrip("/"))
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/vcfannotation/all_samples.{{aligner}}.ann.log"
    shell:
        """
        ABS_PATH=$(pwd)
        
        # Point the config directly to the single genome folder layout
        CONFIG_PATH="$ABS_PATH/{DATABASES_DIR}/snpEff.config"

        snpEff -Xmx{params.ram} -v {params.genome} \
            -c $CONFIG_PATH \
            -dataDir $ABS_PATH/{params.db_dir} \
            -s {output.stats} {input.vcf} | bcftools view -Oz -o {output.vcf}
            
        bcftools index -t {output.vcf}
        """

rule vcf_interactive_query:
    """
    Extracts explicit, flattened tabular data subsets via custom user search metrics.
    """
    input:
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.ann.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{{aligner}}.ann.vcf.gz.tbi"
    output:
        report = f"{READS_DIR}/vcf/query/interactive_report.{{aligner}}.csv"
    params:
        region = lambda w: f"-r {str(config.get('REGION', ''))}" if config.get('REGION') else "",
        include = lambda w: f"-i '{str(config.get('INCLUDE', ''))}'" if config.get('INCLUDE') else "",
        v_type = lambda w: f"-V {str(config.get('V_TYPE', ''))}" if config.get('V_TYPE') else "",
        fmt = "%CHROM,%POS,%REF,%ALT,[%GT],%INFO/ANN\n"
    conda: "../env/dna_vcf.yaml"
    shell:
        """
        mkdir -p {READS_DIR}/vcf/query
        
        echo "# Query Date: $(date)" > {output.report}
        echo "# Parameters: {params.region} {params.include} {params.v_type}" >> {output.report}
        echo "CHROM,POS,REF,ALT,GT,ANNOTATION" >> {output.report}
        
        bcftools query {params.region} {params.include} {params.v_type} \
            -f '{params.fmt}' {input.vcf} >> {output.report}
        """

rule vcf_filter_by_query:
    """
    Generates isolated, sub-selected VCF variant sets using custom criteria splits.
    """
    input:
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.ann.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{{aligner}}.ann.vcf.gz.tbi"
    output:
        vcf = f"{READS_DIR}/vcf/query/{{query_id}}.{{aligner}}.vcf.gz",
        tbi = f"{READS_DIR}/vcf/query/{{query_id}}.{{aligner}}.vcf.gz.tbi"
    params:
        region = lambda w: f"-r {str(config.get('REGION', '')).split('=')[-1]}" if config.get('REGION') else "",
        include = lambda w: f"-i '{str(config.get('INCLUDE', '')).split('=')[-1]}'" if config.get('INCLUDE') else ""
    conda: "../env/dna_vcf.yaml"
    shell:
        """
        bcftools view {params.region} {params.include} {input.vcf} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """