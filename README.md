_Ipomoea hederacea_ bulked segregant analysis workflow

Overview
This repository contains scripts used to conduct a bulked segregant analysis to investigate the genetic basis of extra-lobing in ivyleaf morning glory.

Project structure:
- raw data: reference genome and compressed ORA data
- processed data: all intermediate files (fastqc, trimming, aligning, etc)
- scripts: all scripts and visualizations

Step 1: processing the raw sequencing data
- Done with a .sh file

Step 2: align
- Done with a .sh file

Step 3: bulked segregant analysis, identification of QTLs and candidate genes, GO enrichment, visualization
- The BSA is done in a .Rmd file and submitted as a job using a .sh file

Data availability:
The raw sequence reads supporting the findings of this study have been deposited in the NCBI SRA under BioProject PRJNA1346405.
These data will be publicly released upon publication. 
