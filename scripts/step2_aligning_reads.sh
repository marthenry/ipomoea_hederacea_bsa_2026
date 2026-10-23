#!/bin/bash
#SBATCH --time=14:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --cpus-per-task=15
#SBATCH --job-name="step2_aligning_reads_new_seq_repo"
#SBATCH --mail-user=martin.henry@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --output=step2_aligning_reads.out

#Load modules
module load StdEnv/2023 bwa-mem2/2.2.1 gcc/12.3 samtools/1.20 gatk/4.6.1.0 picard/3.1.0

#Initialize Variables
SEQREPO=../raw_data/reference_genome
REFGENOME=${SEQREPO}/I_hederacea_v1.0.fa
NTHREADS=15

###-----------------------###
###Concatenate fastq files###
###-----------------------###
#Across lanes, but not across paired read direction or bulk
cat ../processed_data/hard_trim/3_lobed_pool_S36_L003_R1_001_val_1.145bp_5prime.fq.gz ../processed_data/hard_trim/3_lobed_pool_S36_L004_R1_001_val_1.145bp_5prime.fq.gz > \
	../processed_data/hard_trim/3_lobed_pool_S36_R1_combined_lanes.fq.gz

cat ../processed_data/hard_trim/3_lobed_pool_S36_L003_R2_001_val_2.145bp_5prime.fq.gz ../processed_data/hard_trim/3_lobed_pool_S36_L004_R2_001_val_2.145bp_5prime.fq.gz > \
	../processed_data/hard_trim/3_lobed_pool_S36_R2_combined_lanes.fq.gz

cat ../processed_data/hard_trim/4_5_lobed_pool_S37_L003_R1_001_val_1.145bp_5prime.fq.gz ../processed_data/hard_trim/4_5_lobed_pool_S37_L004_R1_001_val_1.145bp_5prime.fq.gz > \
	../processed_data/hard_trim/4_5_lobed_pool_S37_R1_combined_lanes.fq.gz

cat ../processed_data/hard_trim/4_5_lobed_pool_S37_L003_R2_001_val_2.145bp_5prime.fq.gz ../processed_data/hard_trim/4_5_lobed_pool_S37_L004_R2_001_val_2.145bp_5prime.fq.gz > \
	../processed_data/hard_trim/4_5_lobed_pool_S37_R2_combined_lanes.fq.gz

#There are 4 files now: one pair for 3-lobed, one pair for 4/5 lobed

##----------------------###
##Index reference genome###
##----------------------###
bwa-mem2 index -p ../processed_data/I_hederacea_v1.0 $SEQREPO/I_hederacea_v1.0.fa

#Align to reference, sort, mark duplicates, and index bam files
while read i
	do
		echo "${i}_R1_combined_lanes.fq.gz"
		bwa-mem2 mem -t $NTHREADS ../processed_data/I_hederacea_v1.0 ../processed_data/hard_trim/${i}_R1_combined_lanes.fq.gz ../processed_data/hard_trim/${i}_R2_combined_lanes.fq.gz > ../processed_data/${i}_align.bwa
		samtools fixmate -@ $NTHREADS -m -O bam ../processed_data/${i}_align.bwa - | samtools sort -@ $NTHREADS - | samtools markdup -@ $NTHREADS -s - ../processed_data/${i}_analysis_ready.bam
		samtools index -b -@ $NTHREADS ../processed_data/${i}_analysis_ready.bam --output ../processed_data/${i}_analysis_ready.bai
	done < ../processed_data/bulk_names.txt
#Duplicates were marked but not removed

###----------------###
###Calling variants###
###----------------###

#Index fasta genome, and create picard dictionary
samtools faidx $REFGENOME
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary -R $REFGENOME -O ${SEQREPO}/I_hederacea_v1.0.dict

# # #Move into processed data
cd ../processed_data

#Iterate through bam files, add sample names, and call variants
for i in $(ls *_analysis_ready.bam)
	do
		echo ${i::-4}
		java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
			I=${i} \
			O=${i::-4}_withsm.bam \
			RGSM=${i::-4} \
			RGPU=unit1 \
			RGPL=ILLUMINA \
			RGLB=${i::-4}
		gatk GetSampleName \
			-I ${i::-4}_withsm.bam \
			-O ${i::-4}_withsm.txt
	done

samtools index -b -@ $NTHREADS 3_lobed_pool_S36_analysis_ready_withsm.bam --output 3_lobed_pool_S36_analysis_ready_withsm.bai

#The 4/5 lobed pool has slightly more sequencing depth, so this is evening out the depth between the pools
samtools view -s .88 4_5_lobed_pool_S37_analysis_ready_withsm.bam -b -o 4_5_lobed_pool_S37_analysis_ready_withsm_88.bam
samtools index -b -@ $NTHREADS 4_5_lobed_pool_S37_analysis_ready_withsm_88.bam --output 4_5_lobed_pool_S37_analysis_ready_withsm_88.bai

for i in $(ls *withsm*.bam)
	do
		prefix=$(cut -d "." -f1 <<< $i)
		gatk HaplotypeCaller -R $REFGENOME \
			--java-options "-Xmx8g" \
			-I ${prefix}.bam \
			-O ${prefix}.g.vcf.gz \
			--native-pair-hmm-threads $NTHREADS \
			-ploidy 2 \
			-ERC GVCF
	done

find -type f -name "*.g.vcf.gz" > gvcf_list.txt

gatk CombineGVCFs \
	-R $REFGENOME \
	--variant 3_lobed_pool_S36_analysis_ready_withsm.g.vcf.gz \
	--variant 4_5_lobed_pool_S37_analysis_ready_withsm_88.g.vcf.gz \
	-O combined.vcf

#Create table for QTLseqR input
gatk VariantsToTable \
   -V combined.vcf \
   -F CHROM -F POS -F REF -F ALT \
   -GF AD -GF DP -GF GQ -GF PL\
   -R $REFGENOME \
   -O qtlseqr_input_subsampled.tsv