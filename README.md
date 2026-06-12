# Bioinformatics_Toolbox
This is a bioinformatics toolbox running on Snakemake platform via conda environment (tested on linux). \
Main purpose: 
1. Partially/fully automate bioinfo analysis using existing packages with less commands
2. Build an expandable scaffolding that can adapt to various job reqirement.

Featuring:
1. Centeralized control via Snakefile.
2. Modularized bioinfo "tools", lower maintenance burdens.
3. Portable environments via snakemake env control.
4. Tolerates single- and pair-end sequences.
4. Neatly organized output folders having everything readily accessible under reads/project_name folder.
5. Project specific config variables, allows switching between projects with a single name change.
6. Stand alone Reference folder, allows ref sharing across projects, saving spaces.
7. Organized logs and project_map, AI agent ready!


### Bookmarks 
1. Setups \
[1.1 Setup Miniconda 3](#1.1-Miniconda-3) \
[1.2 Setup Snakemake](#1.2-Snakemake)
2. Toolbox Functions \
Universal functions \
[2.1 Preparing runinfo.csv](#2.1-Preparing-runinfo.csv) \
[2.2 Download Reference](#2.2-Download-Reference) \
[2.3 QC](#2.3-QC) \
DNA and VCF related functions \
[2.4 DNA Aligner](#2.4-DNA-Aligner) \
[2.5 VCF](#2.5-VCF) \
[2.6 DNA rigidity score](#2.6-DNA-rigidity-score) \
RNAseq related functions \
[2.7 RNA Aligner](#2.7-RNA-aligner) \
[2.8 RNA expression and tools compare](#2.8-RNA-expression) \
[2.9 Expression report analysis Heatmap and function enrichment](#2.9-expression-report-analysis-heatmap-and-function-enrichment) \
Clean up excessive files \
[2.final Clean up](#2.final-Clean-up)
3. Folder Structure \
[3.Folder Structure tree diagram](#3.Folder-Structure-tree-diagram)






# 1. Setups 
<a id="1.1-Miniconda-3"></a>
## 1.1 Miniconda 3
   1. Download latest miniconda 3
   ```bash
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   ```
   2. Install 
   ```bash
      bash Miniconda3-latest-Linux-x86_64.sh
   ```
   3. Reload shell
   ```bash
      source ~/.bashrc
   ```
<a id="1.2-Snakemake"></a>
## 1.2 Snakemake
   1. Download and unpack github package under your working directory.
   2. Make sure environment.yml is in the working dir.
   3. Move terminal to working dir and run:
   ```bash
      conda env create -f snakemake_install.yml -n new_env_name
   ```
   4. Run following command in directory terminal (first time only), to permit execution of run.sh and bcfquery.sh.
   ```bash
      chmod +x run.sh
      chmod +x bcfquery.sh
   ``` 

# 2. Toolbox Functions
You can run:`./ run.sh your_PRJNAME note ` at any time to check available functions.

<a id="2.1-Preparing-runinfo.csv"></a>
## 2.1 Preparing runinfo.csv
   runinfo.csv is dependent in multiple rules to establish the list of sample for analysis.

   ### Online dataset
   Download a runinfo.csv
   1. Create subdirectory: `reads/your_PRJNAME`. 
   2. Copy `reads/Template/suitable_config.yaml` to path above.
   3. Rename it to `config.yaml`
   4. Enter config variable
   ```yaml
            `PRJNUMBER`: "PRJNA257197" # example
            `PRJNAME`: "your_PRJNAME" 
            `DSOURCE`: "SRA" 
            `N`: 10000 # download a N number of reads from origin fastq
   ```
   5. Acitvate snakemake env, e.g. 
   ```bash
   conda activate snakemake
   ```
   6. Run following code in terminal:
   ```bash
   ./run.sh your_PRJNAME runinfo
   ```
   7. sra_runinfo.csv is downloaded to reads/your_PRJNAME/
  
   ### Selecting a subset of online dataset
   1. locate your subset on its website/journal
   2. use terminal command to select the subset by:
   ```bash
   awk 'NR==1 || /^(SRR00[1-5]|SRR007|SRR00[10-11]),/' sra.runinfo.csv > filtered_runinfo.csv 
   ```
   3. rename filter_runinfo.csv to sra_runinfo.csv

   ### Local dataset
   1. Create subdirectory: reads/your_PRJNAME. 
   2. Copy reads/Template/suitable_config.yaml to path above.
   3. Rename it to `config.yaml`
   4. Enter config variable
   ```yaml
            `PRJNUMBER`: "LOCAL" 
            `PRJNAME`: "your_PRJNAME" 
            `DSOURCE`: "LOCAL"
   ```
   5. Place you fastq files under reads/your_PRJNAME/
   6. Acitvate snakemake env, e.g. 
   ```bash
   conda activate snakemake
   ```
   7. Run following code in terminal:
   ```bash
   ./run.sh your_PRJNAME runinfo
   ```
   8. local_runinfo.csv is generated using fastq file names

<a id="2.2-Download-Reference"></a>
## 2.2 Download Reference
  Enter the following config variables:
   ```yaml
   REF: 
      SOURCE: "ncbi" # choose between "ncbi" and "ensembl"
      SPECIES: "zaire_ebolavirus" # latin name for species
      ASSEMBLY: "GCF_000848505.1" # e.g. ncbi: GRCh38, ensembl: R64-2-1
      RELEASE: "2" # e.g. ncbi: p14 (patch #), ensemble: 83
      ACC: "GCF_000848505.1" # leave empty for ensembl
   ```
   
   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME download_ref
  ```
  
   fasta, gff3 files will be downloaded, gtf file will be generated from gff3. \
   File Locations: \
   `refs/SOURCE/SPECIES/ASSEMBLY/RELEASE`



<a id="2.3-QC"></a>
## 2.3 QC
   Run fastQC, mutltiQC and trim_galore on fastq files
   Trim galore requires config variables 
   ```yaml
         `QUALITY`: "30"  # trims quality below 30
         `LENGTH`: "20"   # minimum trimmed sequence length to kept
   ```

   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME qc
  ```

  File Locations: \
   fastQC and multiqQC results: `reads/your_PRJNAME/qc/` \
   trim galore results: `reads/your_PRJNAME/qc_trimmed/`

<a id="2.4-DNA-Aligner"></a>
## 2.4 DNA Aligner
   Align DNA sequences to ref using bwa, bowtie2 or minimap2_sr, then filtered the bam files by mapping quality and SAMflags with SAMtools \
   Can choose from 3 DNA aligner(`bwa`, `bowtie2`, `minimap2_sr`) by config: \
         `ALIGNER`: "bwa" \
   SAMtools view requires config variables 
   ```yaml
         `SAMFLAG`: "3"   
         `MAPQ`: "30"
   ```
   [Picard Samflag Decoder](https://broadinstitute.github.io/picard/explain-flags.html) 

   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME dna_align
  ```
  File Locations: \
   Aligner index:
   `refs/SOURCE/SPECIES/ASSEMBLY/RELEASE/aligner_name/`\
   Bam files:
   `reads/your_PRJNAME/bam/` \
   Filtered bam files:
   `reads/your_PRJNAME/bam/filtered/`

<a id="2.5-VCF"></a>
## 2.5 VCF
   Run VCF calling(freebayes) > create consensus.fa, merge and normalization (bcftools) > build database and annotate (SnpEff) \
   
   Config for freebayes: 
   ```yaml
         `PLOIDY`: "2" # 1 for viral/haploid, 2 for human/diploid)
   ```
   Config for SnpEff 
   ```yaml
         `RAM`: "4g"  # Ram usage for vcf annotation
   ```
   
   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME dna_vcf
  ```

   File locations:
   Raw vcf samples: `reads/your_PRJNAME/vcf/raw/` \
   Consensus.fa files: `reads/your_PRJNAME/vcf/consensus/` \
   Merged and normalized vcf: `reads/your_PRJNAME/vcf/` \
   SnpEff database: `database/snpeff/SPECIES.ASSEMBLY.RELEASE/` \
   Annotated vcf: `reads/your_PRJNAME/vcf/`

   ### Querying Annotated VCF with bcfquery.sh
   this is a shortcut built to run:
   1. vcf_interactive_query, to generate query.csv by bcftools query
   2. vcf_filter_by_query, to generate query.vcf.gz and index by using the same query condition(s)
   #### Quick way to query
   ```bash
    ./bcfquery your_project_name REGION="Chr1:100-500" INCLUDE="QUAL>50" VTYPE="indels"
    # or
    ./bcfquery your_project_name VTYPE="indels"
   ```
   Query condition(s) can be choosen from REGION, INCLUDE and VTYPE (for now)

<a id="2.6-DNA-rigidity-score"></a>
## 2.6 DNA rigidity score
   ### What is it?
   It is a pet project inspired by Gemini AI on the topic of TurboQuant where we found similarity between the this latest spherical compression logic and the mechanism of DNA compression by histone.

   ### How to score DNA rigidity?
   1. Create project folder: `reads/your_PRJNAME` 
   2. Copy dna_rigidity_config.yaml to you project file and rename to config.yaml 
   3. If you are running an ongoing VCF or other DNA project, their config template will cover all for rigidity.

   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME dna_rigid
   ```
   
   File locations:
   tsv and bedgraph: `reads/your_project_name/dna_rigidity/`

   ### Score interpretation: 
   AT rich locale: 5.0 (rigid) \
   GC rich locale: -3.0 (flexible) \
   Baseline/default: 0.0 \
   bedgraph can be piped into IGV for paralleled comparison with other results.


<a id="2.7-RNA-aligner"></a>
## 2.7 RNA Aligner
   Route 1: Bam by sequence alinger (Star, Hisat2) > convert bam to bigwig > FeatureCount > merge \
   Route 2: Pseudo counts by pseudo aligner (Salmon, Kallisto) > merge 

   Can choose from RNA aligners Star, Hisat2, Salmon, Kallisto by config: \
         `RNA_ALIGNER`: "star" \
   Star config variable 
   ```yaml
         `STAR_SA`: "10"   # suffix array index bases, see config instruction for details
   ```
   Kallisto config variable 
   ```yaml
         KALLISTO_FRAG_LEN: "40" # default 200
         KALLISTO_FRAG_SD: "20" # default 20 (%)
   ```
   FeatureCount config variables 
   ```yaml
         `FEATURE_TYPE`: "exon" 
         `ID_ATTR`: "gene_id" # how to group and sum up
   ```

   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME rna_align
  ```
   File Locations: \
   Aligner index:
   `refs/SOURCE/SPECIES/ASSEMBLY/RELEASE/aligner_name/`\ 
   Bam files:
   `reads/your_PRJNAME/bam/` \
   Bigwig files:
   `reads/your_PRJNAME/bam/bigwig/` \
   FeatureCount files:
   `reads/your_PRJNAME/counts/sequence_individual/` \
   Pseudo counts files:
   `reads/your_PRJNAME/counts/pseudo_individual/` \
   Merged counts file:
   `reads/your_PRJNAME/counts/`

<a id="2.8-RNA-expression"></a>
## 2.8 RNA expression and tools compare
   Expression analysis by DESeq2 or edgeR > Compare Stat of the two analysers

   Can choose DESeq2 or edgeR by config: 
   ```yaml
         `EXP_ANALYSER`: "deseq2" \
   ```
 DESeq2 config variable 
   ```yaml
         `DESEQ2_NORM`: "TRUE"   # use DESeq2 normalize method
   ```
   edgeR config variables 
   ```yaml
         `EDGER_NORM`: "TMM"     #TMM, RLE, upperquartile, none
   ```
   
   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME rna_exp
   ```
   Re-run above using a different expression analyser, then use:
   ```bash
   ./run.sh your_PRJNAME tools_compare
   ```
   to get comparison stats and plot. 

   File Locations: \
   expression results: `reads/your_PRJNAME/expression/` \
   tools compare results: `reads/your_PRJNAME/expression/` \
   tools compare plots: `reads/your_PRJNAME/expression/plots/ 

<a id="2.9-Expression-report-analysis-Heatmap-and-function-enrichment"></a>
## 2.9 Expression report analysis, Heatmap and function enrichment
   PCA > MA plot > Volcano plot > P-value histogram > Significant differentially expressed gene CSV > Heatmap > Enrichment dotplot and csv

   Universal config variables: 
   ```yaml
         `RNA_FDR`: 0.05 # adj P-value threshold
         `RNA_LFC`: 0.58 # log fold change threshold
         `BATCH_COLUMN`: "condition" # sample cluster defined by condition column in runinfo.csv
         `BG_COLOR`: "white" # background color for plots: white, transparent/na, any hex code
   ```
   function enrichment dotplot config variable:
   ```yaml
         `DOTPLOT_HEIGHT`: 16   # adjest height of dotplot to resolve label coverlap
   ```
   
   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME rna_report
   `/run.sh your_PRJNAME rna_enrich
   ```

   File Locations: \
   PCA, Volcano, MA plot, Heatmap, P-value histgram, enrichment dotplot: \
   `reads/your_PRJNAME/expression/plots` \
   significant Diff Expression.csv, enrichment GO result: \
    `reads/your_PRJNAME/expression/` 



<a id="2.final-Clean-up"></a>
## 2.final Clean up

   For Snakemake to accept both single and paired end sequences without throwing a "tamtrum", dummy r2 files are created in multiple steps.

   To clean up the dummy files, use command in terminal:
   ```bash
   ./run.sh your_project_name cleanup
   ```
   To re-run `cleanup` function, remove cleanup_complete.done under path: \
   `reads/your_project_name/log/` 


   !!! Only use this function after all analysis were done, snakemake will re-create these files when further analysis was conducted !!!



<a id="3.Folder-Structure-tree-diagram"></a>
## 3.Folder Structure tree diagram

```text
Working Directory 
├── databases
│   └── snpEff                        (snpeff database for vcf annotation)
│       └── {SPECIES}.{ASSEMBLY}.{RELEASE}   (databases built in analysis)
│           ├── sequence.bin
│           ├── snpEff.config
│           └── snpEffectPredictor.bin
├── env
│   └── example_env.yaml              (portable environment configs)
├── reads
│   ├── {RNAseq_project}
│   │   ├── bam                       (bam file and index by Star or Hisat2 aligner)
│   │   │   └── bigwig                (bigwig files for better IGV visualization)
│   │   ├── counts
│   │   │   ├── pseudo_individual     (gene count by kallisto or Salmon pseudo aligner)
│   │   │   └── sequence_individual   (gene count by FeatureCount from bam)
│   │   ├── expression                (expression results from DESeq2, edgeR, fold change and functional enrichment)
│   │   │   └── plots                 (pca, diagnostic plots, volcano, heatmap and enrichment dotplot)
│   │   ├── logs                      (log files for AI agent)
│   │   ├── qc                        (qc and multiqc files)
│   │   ├── qc_trimmed                (trimmed fq.gz files)
│   │   ├── config.yaml               (project specific variables)
│   │   ├── sra/local_runinfo.csv     (sample list)
│   │   └── example_1.fastq
│   │
│   └── {VCF_project}
│       ├── logs                      (log files for AI agent)
│       ├── bam                       (bam file and index)
│       │   └── filtered              (filtered bam files and index)
│       ├── qc                        (fastqc and multiqc files)
│       ├── qc_trimmed                (trimmed fastq files)
│       ├── vcf                       (merged, normalized, annotated vcf files)
│       │   ├── consensus             (consensus fastq files)
│       │   └── query                 (stores query.csv, query.vcf.gz and index)
│       ├── dna_rigidity              (tsv and bedgraph files containing DNA rigidity scores)
│       ├── config.yaml               (project specific variables)
│       ├── sra/local_runinfo.csv     (sample list)
│       ├── sample_1.fastq
│       └── sample_2.fastq
│
├── refs
│     ├── ensembl                     
│     │    └── species
│     │        └── assembly
│     │              └── release      (ensembl reference files, fa, gff3)
│     └── ncbi 
│          └── species
│              └── assembly
│                    └── release      (ncbi reference files, fa, gff3)
│
├── toolbox                           (Modularized tools)
│   └── scripts                       (python and R scripts)
│
├── Snakefile                         (central control)
├── run.sh                            (command shortcut)
├── bcfquery.sh                       (commands to query annotated vcf)
├── PROJECT_MAP.md                    (project map for AI agent)
├── snakemake_install.yml             (portable snakemake env setup file)
└── README.md

```


