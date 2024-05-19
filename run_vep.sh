#!/bin/bash
#20240529
############### preprocessing ###############
# making aliases in ~/.bash_aliases runable in script
shopt -s expand_aliases
source ~/.bash_aliases
# config 
env_name="WES_preprocessing_SNP" # Àô¹Ò¦WºÙ½Ð½T»{
# enable bash script calling conda function
eval "$(conda shell.bash hook)"
conda activate ${env_name}

# ref config
ASSEMBLY_FASTA=/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa
ASSEMBLY=GRCh38

DIR_CACHE=/home/data/database/vep
CACHE_VERSION=103

THREAD=20

# project config
input_vcf="all_samples_2callers.vcf"
output_vcf="all_samples_2callers.vep103.vcf"

############### function ###############
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "start time: $start_time"

vep -i $input_vcf -o $output_vcf --vcf --cache --dir_cache $DIR_CACHE --cache_version $CACHE_VERSION \
--assembly $ASSEMBLY --force_overwrite --fork $THREAD --offline --no_stats --everything -pick \
--fasta $ASSEMBLY_FASTA

end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "end time: $end_time"


