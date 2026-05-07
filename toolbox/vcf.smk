localrules: vcfcall, vcfmerge, vcfnorm, vcfconsensus, download_snpeff_db, vcfannotation

rule vcfcall:
    input:
        bam = f"reads/{PRJNAME}/bam/filtered/{{sample}}.{{aligner}}.filtered.bam",
        ref = f"refs/{config['REFNAME']}/{config['ACC']}.fa",
        fai = f"refs/{config['REFNAME']}/{config['ACC']}.fa.fai"
    output:
        vcf = temp(f"reads/{PRJNAME}/vcf/raw/{{sample}}.{{aligner}}.vcf.gz"),
        tbi = temp(f"reads/{PRJNAME}/vcf/raw/{{sample}}.{{aligner}}.vcf.gz.tbi")
    params:
        ploidy = config.get("PLOIDY", 2)
    conda: "../env/dna_vcf.yaml"
    log: f"reads/{PRJNAME}/logs/vcf/vcfcall/{{sample}}.{{aligner}}.log"
    shell:
        """
        freebayes -f {input.ref} --ploidy {params.ploidy} {input.bam} | bcftools view -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """

rule vcfmerge:
    input:
        vcfs = lambda w: expand(f"reads/{PRJNAME}/vcf/raw/{{sample}}.{w.aligner}.vcf.gz", sample=get_samples(w)),
        tbis = lambda w: expand(f"reads/{PRJNAME}/vcf/raw/{{sample}}.{w.aligner}.vcf.gz.tbi", sample=get_samples(w))
        #create a list of filepaths e.g. (reads/PRJNAME/vcf/raw/SRR001, reads/PRJNAME/vcf/raw/SRR002)
    output:
        merged_vcf = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.merged.vcf.gz"
    conda: "../env/dna_vcf.yaml"
    log: f"reads/{PRJNAME}/logs/vcf/vcfmerge/all_sample.{{aligner}}.merged.log"
    shell:
        "bcftools merge {input.vcfs} -Oz -o {output.merged_vcf}"

rule vcfnorm:
    input:
        vcf = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.merged.vcf.gz",
        ref = f"refs/{config['REFNAME']}/{config['ACC']}.fa"
    output:
        vcf_gz = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz",
        tbi = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz.tbi"
    conda: "../env/dna_vcf.yaml"
    log: f"reads/{PRJNAME}/logs/vcf/vcfnorm/all_samples.{{aligner}}.merged.norm.log"
    shell:
        """
        bcftools norm -f {input.ref} -m -any {input.vcf} -Oz -o {output.vcf_gz}
        bcftools index -t {output.vcf_gz}
        """

rule vcfconsensus:
    input:
        vcf = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz",
        tbi = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz.tbi",
        ref = f"refs/{config['REFNAME']}/{config['ACC']}.fa"
    output:
        fasta = f"reads/{PRJNAME}/vcf/consensus/{{sample}}.{{aligner}}.consensus.fa"
    conda: "../env/dna_vcf.yaml"
    log: f"reads/{PRJNAME}/logs/vcf/vcfconsensus/{{sample}}.{{aligner}}.consensus.log"
    shell:
        "bcftools consensus -f {input.ref} -s {wildcards.sample} {input.vcf} > {output.fasta}"

# this rule was placed on hold because building a local db is more desirable 
# rule download_snpeff_db:
#     output:
#         directory(f"databases/snpeff/{config['SNPEFF_DB']}")
#     params:
#         db = config.get("SNPEFF_DB", "EBOLA_VIRUS_GENOME"),
#         ram = config.get("RAM", "4g")
#     conda: "../env/dna_vcf.yaml"
#     log: f"reads/{PRJNAME}/logs/vcf/snpeff_download.log"
#     shell:
#         """
#         snpEff -Xmx{params.ram} download -v -dataDir databases/snpeff {params.db} > {log} 2>&1
#         """

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
    log: f"reads/{PRJNAME}/logs/vcf/snpeff_build_{config['REFNAME']}.log"
    shell:
        """
        # Step 1: create subdirectory, copy fa and gff3 files
        mkdir -p {params.db_dir}/{params.genome}
        cp {input.fasta} {params.db_dir}/{params.genome}/sequences.fa
        cp {input.gff} {params.db_dir}/{params.genome}/genes.gff

        # Step 2: create local snpeff.config
        # pwd is used to prevent SnpEff from creating "double path"
        # e.g. SnpEff will create databases/snfeff/refname/databases/snpeff/refname
        # because it recognize snpeff.config location as the working direct (which is under databases/snfeff/refname)
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
        vcf = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.merged.norm.vcf.gz",
        db = f"databases/snpeff/{config['REFNAME']}/snpEffectPredictor.bin" 
    output:
        vcf = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.ann.vcf.gz",
        tbi = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.ann.vcf.gz.tbi",
        stats = f"reads/{PRJNAME}/vcf/all_samples.{{aligner}}.snpeff_stats.html"
    params:
        genome = config['REFNAME'],
        ram = config.get("RAM", "4g"),
        config_path = f"databases/snpeff/{config['REFNAME']}/snpEff.config",
        db_dir = "databases/snpeff"
    conda: "../env/dna_vcf.yaml"
    log: f"reads/{PRJNAME}/logs/vcf/vcfannotation/all_samples.{{aligner}}.ann.log"
    shell:
        """
        # Capture current directory to prevent path doubling
        ABS_PATH=$(pwd)
        snpEff -Xmx{params.ram} ann \
            -c $ABS_PATH/{params.config_path} \
            -dataDir $ABS_PATH/{params.db_dir} \
            -s {output.stats} \
            {params.genome} \
            {input.vcf} | bcftools view -Oz -o {output.vcf}
        
        # Index the new file to create the .tbi
        bcftools index -t {output.vcf}
        """
        
rule vcf_interactive_query:
    input:
        vcf = f"reads/{PRJNAME}/vcf/all_samples.{config['ALIGNER']}.ann.vcf.gz",
        tbi = f"reads/{PRJNAME}/vcf/all_samples.{config['ALIGNER']}.ann.vcf.gz.tbi"
    output:
        report = f"reads/{PRJNAME}/vcf/query/{{query_id}}.csv"
    params:
        region = lambda w: f"-r {str(config.get('REGION', '')).split('=')[-1]}" if config.get('REGION') else "",
        include = lambda w: f"-i '{str(config.get('INCLUDE', '')).split('=')[-1]}'" if config.get('INCLUDE') else "",
        v_type = lambda w: f"-i 'TYPE=\"{str(config.get('VTYPE', '')).split('=')[-1]}\"'" if config.get('VTYPE') else "",

        fmt = config.get("FMT", "%CHROM,%POS,%REF,%ALT,[%GT],%INFO/ANN\\n")
    conda: "../env/dna_vcf.yaml"
    shell:
        """
        mkdir -p reads/{PRJNAME}/vcf/query
        
        # Header with metadata for future reference
        echo "# Query Date: $(date)" > {output.report}
        echo "# Parameters: {params.region} {params.include} {params.v_type}" >> {output.report}
        echo "CHROM,POS,REF,ALT,GT,ANNOTATION" >> {output.report}
        
        bcftools query {params.region} {params.include} {params.v_type} \
            -f '{params.fmt}' {input.vcf} >> {output.report}
        """

rule vcf_filter_by_query:
    input:
        vcf = f"reads/{PRJNAME}/vcf/all_samples.{config['ALIGNER']}.ann.vcf.gz",
        tbi = f"reads/{PRJNAME}/vcf/all_samples.{config['ALIGNER']}.ann.vcf.gz.tbi"
    output:
        vcf = f"reads/{PRJNAME}/vcf/query/{{query_id}}.vcf.gz",
        tbi = f"reads/{PRJNAME}/vcf/query/{{query_id}}.vcf.gz.tbi"
    params:
        region = lambda w: f"-r {str(config.get('REGION', '')).split('=')[-1]}" if config.get('REGION') else "",
        include = lambda w: f"-i '{str(config.get('INCLUDE', '')).split('=')[-1]}'" if config.get('INCLUDE') else "",
        # bcftools view uses -v indels directly
        v_type = lambda w: f"-v {str(config.get('VTYPE', '')).split('=')[-1]}" if config.get('VTYPE') else ""
    conda: "../env/dna_vcf.yaml"
    shell:
        """
        bcftools view {params.region} {params.include} {params.v_type} -O z -o {output.vcf} {input.vcf}
        bcftools index -t {output.vcf}
        """