#!/bin/bash
#20240529
############### preprocessing ###############
# config 
env_name="WES_preprocessing_SNP" # 環境名稱請確認
# enable bash script calling conda function
eval "$(conda shell.bash hook)"
conda activate ${env_name}

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet_test.2024-02-26.tsv"

# Create output and log directories if they don't exist
input_fastq_dir="./fastq"
output_fastq_dir="./fastq"

############### function ###############
run_change_id(){
    file_id=${1}
    Sample_id=${2}
    mv ${input_fastq_dir}/${file_id}.fastq.gz ${output_fastq_dir}/${Sample_id}.fastq.gz
}

############### main ###############
tail -n +2 $sample_sheet | while IFS=$'\t' read -r -a sample
do
    file_ID=${sample[0]} # [0]: file ID
    TCGA_ID=${sample[6]} # [6]: sample ID
    run_change_id ${file_ID} ${TCGA_ID}
done



