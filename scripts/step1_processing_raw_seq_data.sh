#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --cpus-per-task=8
#SBATCH --job-name="step1_processing_raw_seq_data"
#SBATCH --mail-user=martin.henry@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --output=step1_processing_raw_seq_data.out

#Load environments
module load StdEnv/2023 fastqc/0.12.1 python/3.13.2 arrow/22.0.0

#Initialize Variables
SEQREPO=~/projects/def-jstinchc/sequencing_repo/Ipomoea_hederacea/Ihederacea_five_lobed_2025 #location of raw sequencing data in ora format
ORADPATH=~/projects/def-jstinchc/code_repo/orad.2.7.0.linux/orad #path to orad executable
ORAREF=../raw_data/oradata/homo_sapiens #location of ora fastq reference
FASTQREPO=../processed_data/fastq_raw #where to store decompressed fastq files
# MULTIQC=~/multiqc_env/bin/activate #location of multiqc environment
TRIMGALORE=~/TrimGalore/bin/activate #location of TrimGalore environment
trim_galore=~/projects/def-jstinchc/mhenry29/TrimGalore-0.6.10/trim_galore #executable for TrimGalore
NTHREADS=8



# #--------------------------#
# ###Decompressing raw data###
# #--------------------------#

# #Create a folder where you will store the .ora files that are the original raw sequencing data
# #mkdir $SEQREPO
# #Create a folder where the fastq files go
# mkdir $FASTQREPO

# #Download and decompress orad references
# #wget -O ../raw_data/ora_allother.refbin.tar.gz https://webdata.illumina.com/downloads/software/dragen/references/ora-270/oradata_all_species.tar.gz
# #tar -xzf ../raw_data/ora_allother.refbin.tar.gz -C ../raw_data/

# #Decompress the data
# for i in $(ls ${SEQREPO}/*.ora) 
# do
#     $ORADPATH --gz -P $FASTQREPO --ora-reference $ORAREF $i
# done

# #Send output to results
# mkdir ../processed_data/fastqc_pre_trim

# #Do the fastqc
# fastqc -o ../processed_data/fastqc_pre_trim ${FASTQREPO}/*.gz

# #Activate virtual environment to do multiqc next
# source $MULTIQC

# #Run multiqc
# mkdir ../processed_data/multiqc_pre_trim
# multiqc ../processed_data/fastqc_pre_trim --outdir ../processed_data/multiqc_pre_trim

# #-----------------------#
# ###Soft trimming reads###
# #-----------------------#

#Trimming FastQ files using TrimGalore
##Go into virtual environment -> TrimGalore
source $TRIMGALORE

#Install cutadapt in virtual environment
#pip install cutadapt

#Run trim_galore variable and do a while loop for paired-end files
#20 is used for --2colour because it is the default phred score for the -q, which is the standard quality cutoff
while read i
do
echo "start trimming $i"
$trim_galore --fastqc --2colour 20 --cores $NTHREADS -o ../processed_data/soft_trim --paired ${FASTQREPO}/$i\_R1_001.fastq.gz ${FASTQREPO}/$i\_R2_001.fastq.gz
echo "done trimming $i"
done<../processed_data/read_names.txt

multiqc ../processed_data/soft_trim --outdir ../processed_data/soft_trim/trimmed_multiqc_output

#-----------------------#
###Hard trimming reads###
#-----------------------#

#Run trim_galore variable and do a while loop for paired-end files
#20 is used for --2colour because it is the default phred score for the -q, which is the standard quality cutoff
while read i
do
echo "start trimming $i"
$trim_galore --fastqc --hardtrim5 145 --cores $NTHREADS -o ../processed_data/hard_trim \
--paired ../processed_data/soft_trim/$i\_R1_001_val_1.fq.gz \
../processed_data/soft_trim/$i\_R2_001_val_2.fq.gz
echo "done trimming $i"
done<../processed_data/read_names.txt

multiqc ../processed_data/hard_trim --outdir ../processed_data/hard_trim/trimmed_multiqc_output

