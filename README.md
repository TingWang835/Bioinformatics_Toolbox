# Bioinfomatics-Project
## Snakemake bioinfomatics toolbox
This is a bioinfomatics toolbox running on Snakemake platform via conda environment (tested on linux).
Main purposse: to partially/fully automate bioinfo analysis using existing packages but with less commands

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
4. run following command in terminal (first time only), to allow excution of run.sh and bcfquery.sh
```bash
 chmod +x run.sh
 chmod +x bcfquery.sh
``` 

## How to run VCF
### SRA online project
   To analyse data requires downloading from SRA: 

   1. create subdirectory: reads/your_project_name. 
   2. copy reads/PRJNAME/config.yaml to path above.
   3. enter your project specific variables, make sure:
   ```bash
    DSOURCE: "SRA"
   ```
   4. activate your snakemake environment e.g.   
   ```bash
    conda activate snakemake
   ```  
   5. run the following command in terminal 
   ```bash
   ./run.sh your_project_name note 
   ```
   6. select from functions listed, replacing "note" e.g.
   ```bash
   ./run.sh your_project_name vcf_all 
   ```

### Local project
   To analyse data from your local sequencer: 

   1. create subdirectory: reads/your_project_name. 
   2. copy local fq/fastq/fq.gz file(s) to subdirectory above, a local_runinfo.csv will be generated based on sample names.
   3. copy reads/PRJNAME/config.yaml to the same path.
   4. enter your project specific variables, make sure:
   ```bash
    DSOURCE: "LOCAL"
   ``` 
   5. activate your snakemake environment e.g.   
   ```bash
    conda activate snakemake
   ```  
   6. run the following command in terminal 
   ```bash
   ./run.sh your_project_name note 
   ```
   7. select from functions listed, replacing "note" e.g.
   ```bash
   ./run.sh your_project_name vcf_all 
   ```

### Running selected samples
  if you would like to analysis a selected group(s) of samples, remove the unwanted rows from local/sra_runinfo.csv, sample name(s) listed in runinfo.csv will determine which sample(s) will be analysed.


### bcfquery.sh
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
│   ├── aligner.smk                (align trimmed fastq to refs to create bam and filtered bam)
│   ├── getdata.smk                (download/register fastq under reads/PRJNAME, download fa and gff from SRA and create index)
│   ├── qc.smk                     (QC fastq files and trim)
│   └── vcf.smk                    (call, merge, norm and annotate vcf, also creates consensus fastq)
│
├── Snakefile                      (central control)
├── run.sh                         (command shortcut)
├── bcfquery.sh                    (commands to query annotated vcf)
├── PROJECT_MAP.md                 (project map for AI agent)
└── README.md

```


