#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=50000M
#SBATCH --cpus-per-task=1
#SBATCH --job-name="step3_bulked_segregant_analysis"
#SBATCH --mail-user=martin.henry@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --output=step3_bulked_segregant_analysis.out

module load StdEnv/2023 r/4.4.0

Rscript -e "rmarkdown::render(input = 'step3_bulked_segregant_analysis.Rmd')"
