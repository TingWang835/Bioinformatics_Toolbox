localrules: vcfcall, vcfmerge, vcfnorm, vcfconsensus, snpeff_build, vcfannotation, vcf_interactive_query, vcf_filter_by_query

# =============================================================================
# Rules
# =============================================================================

rule vcfcall:
    input:
        unpack(get_ref_source),
        bam = f"{READS_DIR}/bam/filtered/{{sample}}.{{aligner}}.filtered.bam",
        fai = lambda wildcards: f"{get_ref_source(wildcards)['ref']}.fai"
    output:
        vcf = temp(f"{READS_DIR}/vcf/raw/{{sample}}.{{aligner}}.vcf.gz"),
        tbi = temp(f"{READS_DIR}/vcf/raw/{{sample}}.{{aligner}}.vcf.gz.tbi")
    params:
        ploidy = config.get("PLOIDY", 2)
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/vcfcall/{{sample}}.{{aligner}}.log"
    shell:
        """
        freebayes -f {input.ref} --ploidy {params.ploidy} {input.bam} | bcftools view -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """

rule vcfmerge:
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
    input:
        unpack(get_ref_source),
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.vcf.gz"
    output:
        vcf_gz = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz.tbi"
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/vcfnorm/all_samples.{{aligner}}.merged.norm.log"
    shell:
        """
        bcftools norm -f {input.ref} -m -any {input.vcf} -Oz -o {output.vcf_gz}
        bcftools index -t {output.vcf_gz}
        """

rule vcfconsensus:
    input:
        unpack(get_ref_source),
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz.tbi"
    output:
        fasta = f"{READS_DIR}/vcf/consensus/{{sample}}.{{aligner}}.consensus.fa"
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/vcfconsensus/{{sample}}.{{aligner}}.consensus.log"
    shell:
        "bcftools consensus -f {input.ref} -s {wildcards.sample} {input.vcf} > {output.fasta}"

# =============================================================================
# Functional Annotation & Queries (Pending Refactoring)
# =============================================================================

rule snpeff_build:
    input:
        fasta = f"refs/{config['REFNAME']}/{config['ACC']}.fa",
        gff = f"refs/{config['REFNAME']}/{config['ACC']}.gff"
    output:
        db = f"databases/snpeff/{config['REFNAME']}/snpEffectPredictor.bin"
    params:
        genome = config['REFNAME'],
        db_dir = "databases/snpeff"
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/snpeff_build_{config['REFNAME']}.log"
    shell:
        """
        # Step 1: create subdirectory, copy fa and gff3 files
        mkdir -p {params.db_dir}/{params.genome}
        cp {input.fasta} {params.db_dir}/{params.genome}/sequences.fa
        cp {input.gff} {params.db_dir}/{params.genome}/genes.gff

        # Step 2: create local snpeff.config
        ABS_PATH=$(pwd)
        echo "{params.genome}.genome : {params.genome}" > {params.db_dir}/{params.genome}/snpEff.config

        # Step 3: building DB
        snpEff build -gff3 -v {params.genome} \
            -c $ABS_PATH/{params.db_dir}/{params.genome}/snpEff.config \
            -dataDir $ABS_PATH/{params.db_dir} \
            -noCheckCds -noCheckProtein -nodownload > {log} 2>&1

        # Step 4: cleaning up
        if [ -f {output.db} ]; then
            rm {params.db_dir}/{params.genome}/sequences.fa
            rm {params.db_dir}/{params.genome}/genes.gff
        else
            exit 1
        fi
        """

rule vcfannotation:
    input:
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz",
        db = f"databases/snpeff/{config['REFNAME']}/snpEffectPredictor.bin"
    output:
        vcf = f"{READS_DIR}/vcf/all_samples.{{aligner}}.ann.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{{aligner}}.ann.vcf.gz.tbi",
        stats = f"{READS_DIR}/vcf/all_samples.{{aligner}}.snpeff_stats.html"
    params:
        genome = config['REFNAME'],
        ram = config.get("RAM", "4g"),
        db_dir = "databases/snpeff"
    conda: "../env/dna_vcf.yaml"
    log: f"{LOG_DIR}/vcf/vcfannotation/all_samples.{{aligner}}.ann.log"
    shell:
        """
        ABS_PATH=$(pwd)
        snpEff -Xmx{params.ram} -v {params.genome} \
            -c $ABS_PATH/{params.db_dir}/{params.genome}/snpEff.config \
            -dataDir $ABS_PATH/{params.db_dir} \
            -s {output.stats} {input.vcf} | bcftools view -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """

rule vcf_interactive_query:
    input:
        vcf = f"{READS_DIR}/vcf/all_samples.{dna_aligner}.ann.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{dna_aligner}.ann.vcf.gz.tbi"
    output:
        report = f"{READS_DIR}/vcf/query/interactive_report.csv"
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
    input:
        vcf = f"{READS_DIR}/vcf/all_samples.{dna_aligner}.ann.vcf.gz",
        tbi = f"{READS_DIR}/vcf/all_samples.{dna_aligner}.ann.vcf.gz.tbi"
    output:
        vcf = f"{READS_DIR}/vcf/query/{{query_id}}.vcf.gz",
        tbi = f"{READS_DIR}/vcf/query/{{query_id}}.vcf.gz.tbi"
    params:
        region = lambda w: f"-r {str(config.get('REGION', '')).split('=')[-1]}" if config.get('REGION') else "",
        include = lambda w: f"-i '{str(config.get('INCLUDE', '')).split('=')[-1]}'" if config.get('INCLUDE') else ""
    conda: "../env/dna_vcf.yaml"
    shell:
        """
        bcftools view {params.region} {params.include} {input.vcf} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """