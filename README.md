# Bioinfomatics-Project
## Snakemake bioinfomatics toolbox
This is a bioinfomatics toolbox running on Snakemake platform via conda environment (tesed on linux).

Featuring:
1. Centeralized control via snakefile
2. Modularized bioinfo "tools", runs specific task, easy to maintain
3. Portable and readily mangable environments
4. Stand alone reference folder, allows ref sharing across projects
5. Project specific config.yaml, records variables, allows quick and easy switch between projects
6. Organized logs and project_map, AI agent ready!


## How to Setup
1. Download and unpack github package under your working directory 
2. install miniconda 3 
3. install snakemake via miniconda 3

## How to run
### SRA online project
   To analyse data requires downloading from SRA: 

   1. create subdirectory: reads/your_project_name. 
   2. copy reads/PRJNAME/config.yaml to path above.
   3. enter your project specific variables. 
   4. activate your snakemake environment e.g.   
   ```bash
    conda activate snakemake
   ```  
   5. run the following command in terminal 
   ```bash
   bash run.sh your_project_name note 
   ```
   6. select the function listed replacing "note" e.g.
   ```bash
   bash run.sh your_project_name vcf_all 
   ```

### local project
   To analyse data from your local sequencer: 

   1. create subdirectory: reads/your_project_name. 
   2. copy local fq/fastq/fq.gz file(s) to subdirectory above, a local_runinfo.csv will be generated based on sample names.
   3. copy reads/PRJNAME/config.yaml to the same path.
   4. enter your project specific variables. 
   5. activate your snakemake environment e.g.   
   ```bash
    conda activate snakemake
   ```  
   6. run the following command in terminal 
   ```bash
   bash run.sh your_project_name note 
   ```
   7. select the function listed replacing "note" e.g.
   ```bash
   bash run.sh your_project_name vcf_all 
   ```

### Running selected samples
  if you would like to analysis a selected group(s) of samples, delete the unwanted rows from local/sra_runinfo.csv, sample name(s) listed in runinfo.csv will determine which sample(s) will be analysed.

### VCF query
   


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
8. SAMFLAG: sam flag number for bam filtering, default = 3 (paired + proper paired), decoding by [picard samflag](https://broadinstitute.github.io/picard/explain-flags.html)
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
│   └── example.yaml
│ 
├── reads
│   └── [PRJNAME]
│       ├── logs                   (log files for AI agent)
│       ├── bam                    (bam file and index)
│       │   └── filtered           (filtered bam files and index)
│       ├── qc                     (fastqc and multiqc files)
│       ├── qc_trimmed             (trimmed fastq files)
│       ├── vcf                    (merged, normalized, annotated vcf files)
│       │   └── consensus          (consensus fastq files)
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
├── toolbox                        (modularized tools)
│   ├── aligner.smk
│   ├── getdata.smk
│   ├── qc.smk
│   └── vcf.smk
│
├── Snakefile                      (central control)
├── run.sh                         (command shortcut)
├── bcfquery.sh                    (commands to query annotated vcf)
├── PROJECT_MAP.md                 (project map for AI agent)
└── README.md

```


