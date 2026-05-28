# Bioinformatics_Toolbox
This is a bioinformatics toolbox running on Snakemake platform via conda environment (tested on linux).
Main purposse: 
1. Partially/fully automate bioinfo analysis using existing packages with less commands
2. Build an expandable scaffolding that can adapt to various job reqirement.

Featuring:
1. Centeralized control via snakefile.
2. Modularized bioinfo "tools", lower maintenance burdens.
3. Portable environments via snakemake env control.
4. Neatly organized output folders having everything readily assessible under reads/project_name folder.
5. Project specific config variables, allows switching between projects with a single name change.
6. Stand alone Reference folder, allows ref sharing across projects, saving spaces.
7. Organized logs and project_map, AI agent ready!

Bookmarks
[Chapter 2.2 Download Reference](#2.2-Download-Reference)

# 1. How to Setup 
## 1.1 Setup Miniconda 3
   1. Downlaod latest miniconda 3
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

## 1.2 Setup Snakemake
   1. Download and unpack github package under your working directory.
   2. Makesure environment.yml is in the working dir.
   3. Move terminal to working dir and run:
   ```bash
      conda env create -f snakemake_install.yml -n new_env_name
   ```
   4. Run following command in directory terminal (first time only), to permit excution of run.sh and bcfquery.sh.
   ```bash
      chmod +x run.sh
      chmod +x bcfquery.sh
   ``` 

# 2. Toolbox Functions
## 2.1 Preparing runinfo.csv
   runinfo.csv is dependent in multiple rules to establish the list of sample for analysis.

   ### Online dataset
   Download a runinfo.csv
   1. Create subdirectory: reads/your_PRJNAME. 
   2. Copy reads/Template/suitable_config.yaml to path above.
   3. Rename it to `config.yaml`
   4. Enter \
            `PRJNUMBER`: "PRJNA257197" (example)\
            `PRJNAME`: "your_PRJNAME" \
            `DSOURCE`: "SRA"
   5. Acitvate snakemake env, e.g. 
   ```bash
   conda activate snakemake
   ```
   6. Run following code in terminal:
   ```bash
   ./run.sh your_PRJNAME runinfo
   ```
   7. sra_runinfo.csv is downloaded to reads/your_PRJNAME/
  
   ### Sellecting a subset of online dataset
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
   4. Enter \
            `PRJNUMBER`: "LOCAL" \
            `PRJNAME`: "your_PRJNAME" \
            `DSOURCE`: "LOCAL"
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

<a id="2.2 Download Reference"></a>
## 2.2 Download Reference
   NCBI and Ensembl reference source can be selected.
   
   ### NCBI
   Enter config variables \
          `REFNAME`: "your_REFNAME" \
         `REF_SOURCE`: "NCBI" \
         `ACC`: e.g. "AF086833" (example) 

   ### Ensembl
   Enter config variables \
         `REFNAME`: "your_REFNAME" \
         `REF_SOURCE`: "ENSEMBL" \
         `SPECIES_LATIN`: "saccharomyces_cerevisiae" (example) \
         `ASSEMBLY`: "R64-1-1" (example) \
         `RELEASE`: 112 (example)

   
   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME download_ref
  ```



## 2.2 Running VCF projects
### SRA online project
   To analyse data requires downloading from SRA: 

   1. Create subdirectory: reads/your_project_name. 
   2. Copy reads/PRJNAME/config.yaml to path above.
   3. Enter your project specific variables, make sure:
   ```bash
    DSOURCE: "SRA"
   ```
   4. Activate your snakemake environment e.g.   
   ```bash
    conda activate snakemake
   ```  
   5. Run the following command in terminal. 
   ```bash
   ./run.sh your_project_name note 
   ```
   6. Select from functions listed, replacing "note" e.g.
   ```bash
   ./run.sh your_project_name vcf_all 
   ```

### Local project
   To analyse data from your local sequencer: 

   1. create subdirectory: reads/your_project_name. 
   2. copy local fastq file(s) (accepting fq/fq.gz formats) to subdirectory above.
   3. copy reads/PRJNAME/config.yaml to the same path.
   4. enter your project specific variables, make sure:
   ```bash
    DSOURCE: "LOCAL"
   ``` 
   5. activate your snakemake environment e.g.   
   ```bash
    conda activate snakemake
   ```  
   6. run the following command in terminal.
   ```bash
   ./run.sh your_project_name note 
   ```
   7. select from functions listed, replacing "note" e.g.
   ```bash
   ./run.sh your_project_name vcf_all 
   ```

### Running selected samples
  This toolbox uses local/sra_runinfo.csv as a list to determine which sample will be analysed.

  In a SRA project, runinfo.csv will be downloaded using the PRJNMUBER variable.

  Whereas in local project, a runinfo.csv will be generated based on the sample name(s) of the fastq files under reads/project_name folder.

  To run an analysis on a selected group(s) of samples, remove the unwanted rows from local/sra_runinfo.csv, and run the toolbox as shown above.


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
   Query condition can be choosen from REGION, INCLUDE and VTYPE (for now)

### Cleanup
   For Snakemake to accept both single and paired end sequences without throwing a "tamturm", dummy r2 files are created in multiple steps.

   To clean up the dummy files, use command in terminal:
   ```bash
   ./run.sh your_project_name cleanup
   ```
   !!! Only use this function after all analysis were done, snakemake will re-create these files when further analysis was conducted !!!



## 2.2 DNA rigidity score
### What is it?
It is a pet project inspired by Gemini 3.0 when we chatted about TurboQuant where we found similarity between the this latest spherical compression logic and the mechanism of DNA compression by histone.

### How to score DNA rigidity?
   1. Create subdirectory: reads/your_project_name. 
   2. Copy reads/PRJNAME/config.yaml to path above.
   3. Enter the variables for downloading the right reference. (ACC and REFNAME)
   4. activate your snakemake environment e.g.   
   ```bash
    conda activate snakemake
   ```  
   5. run the following command in terminal.
   ```bash
   ./run.sh your_project_name rigid 
   ```
   6. tsv and bedgraph files will be produced with chromosome as id, position and rigidity score. For AT rich rigid position the score is 5.0, and GC rich flexible position: -3.0, baseline/default: 0.0, bedgraph can be piped into IGV for paralleled comparison with other results.






# 3. Variables in config.yaml (What they do)
### Common variables
1. ACC: identify reference (fa and gff3) for download; 
2. REFNAME: naming subdirectory for reference files under refs/
3. PRJNUMBER: Identify online project from SRA, for accquiring runinfo.csv
4. PRJNAME: naming subdirectory for sample(s) and other analysis outputs under reads/

### getdata.smk specific variables
5. DSOURCE: data source, select between 'SRA' and 'LOCAL'
6. N: the number of reads grep from online SRR.fastq files, default = 10000

### aligner.smk specific variables
7. ALIGNER: select aligner package (available: bwa, bowtie2, minimap2_sr)
8. SAMFLAG: sam flag number for bam filtering, default = 3 (paired + proper paired), encoding by [picard samflag](https://broadinstitute.github.io/picard/explain-flags.html)
9. MAPQ: mapping quality score for bam filtering

### vcf.smk specific variables
10. PLOIDY: viral/haploid: 1, human/diploid: 2, default = 2
11. RAM: ram space used for vcf annotation, default = 4g



## Folder Structure
```text
Working Directory 
├── databases
│   └── snpEff                        (snpeff database for vcf annotation)
│       └── {REFNAME}                 (databases built in analysis)
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
│   └── {REFNAME}
│       ├── ensembl                   (ncbi reference files, fa, gff3)
│       └── ncbi                      (ncbi reference files, fa, gff3)
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


