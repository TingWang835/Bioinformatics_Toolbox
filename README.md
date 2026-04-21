# Bioinformatics_Toolbox
This is a bioinformatics toolbox running on Snakemake platform via conda environment (tested on linux).
Main purposse: 
1. Partially/fully automate bioinfo analysis using existing packages but with less commands
2. Build an expandable scaffolding that can adapt to various job reqirement.

Featuring:
1. Centeralized control via snakefile.
2. Modularized and specialised bioinfo "tools", lower maintenance burdens.
3. Portable and mangable environments via snakemake env control.
4. Neatly organized output folders having everything readily assessible under reads/project_name folder.
5. Project specific config.yaml, records variables, allows quick switch between projects.
6. Stand alone Reference folder, allows ref sharing across projects, saving spaces, easier management.
7. Organized logs and project_map, AI agent ready!


## Setup Miniconda 3
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

## Setup Snakemake
   1. Download and unpack github package under your working directory.
   2. Makesure environment.yml is in the working dir.
   3. Move terminal to working dir and run:
   ```bash
      conda env create -f environment.yml -n new_env_name
   ```
   4. Run following command in directory terminal (first time only), to permit excution of run.sh and bcfquery.sh.
   ```bash
      chmod +x run.sh
      chmod +x bcfquery.sh
   ``` 

## How to run VCF in the toolbox
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
#### How to run it
```bash
  ./bcfquery your_project_name REGION="Chr1:100-500" INCLUDE="QUAL>50" VTYPE="indels"
  # or
  ./bcfquery your_project_name VTYPE="indels"
 ```
   Query condition can be choosen from REGION, INCLUDE and VTYPE (for now)



## Variables in config.yaml (What they do)
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
├── databases                     (databases built in analysis)
│   └── snpeff
│       └── [REFNAME]
│           ├── sequence.bin
│           ├── snpEff.config
│           └── snpEffectPredictor.bin
│
├── env                            (portable environment configs)
│   └── example_env.yaml
│ 
├── reads
│   └── [PRJNAME]
│       ├── logs                   (log files for AI agent)
│       ├── bam                    (bam file and index)
│       │   └── filtered           (filtered bam files and index)
│       ├── qc                     (fastqc and multiqc files)
│       ├── qc_trimmed             (trimmed fastq files)
│       ├── vcf                    (merged, normalized, annotated vcf files)
│       │   ├── consensus          (consensus fastq files)
│       │   └── query              (stores query.csv, query.vcf.gz and index)
│       ├── config.yaml            (project specific variables)
│       ├── sra/local_runinfo.csv  (sample list)
│       ├── sample_1.fastq
│       └── sample_2.fastq
│
├── refs
│   └── [REFNAME]
│       ├── ACC.fa
│       ├── ACC.fa.fai
│       ├── ACC.gff
│       └── aligner.index
│
├── toolbox                        (Store modularized tools)
│   ├── getdata.smk                (download/register fastq under reads/PRJNAME, download fa and gff from SRA and create index)
│   ├── qc.smk                     (QC fastq files and trim)
│   ├── aligner.smk                (align trimmed fastq to refs to create bam and filtered bam)
│   └── vcf.smk                    (call, merge, norm and annotate vcf, also creates consensus fastq)
│
├── Snakefile                      (central control)
├── run.sh                         (command shortcut)
├── bcfquery.sh                    (commands to query annotated vcf)
├── PROJECT_MAP.md                 (project map for AI agent)
├── environment.yml                (used in snakemake setup)
└── README.md

```


