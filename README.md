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


## Folder Structure
```text
Working Directory 
├── databases                     (databases built in analysis)
│   └── snpeff
│       └── {REFNAME}
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

## How to Setup
1. Download and unpack github package under your working directory 
2. install miniconda 3 
3. install snakemake via miniconda 3

## How to run
### SRA online project
   To analysis data requires downloading from SRA: \
   a. create subdirectory: reads/[your_project_name] \
   b. copy reads/[PRJNAME]/config.yaml to path above\
   c. enter your project specific variables \
   d. activate your snakemake environment e.g.:   
    ```bash
    conda activate snakemake
    ```  
   e. run in terminal: 
   ```bash
   bash run.sh [your_project_name] vcf_all
   ```
