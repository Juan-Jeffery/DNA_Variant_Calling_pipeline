#!/bin/bash
#20240529
############### preprocessing ############### 
# config 
env_name="WES_preprocessing_SNP" # 環境名稱請確認
# enable bash script calling conda function
eval "$(conda shell.bash hook)"
conda activate ${env_name}
threads=20

# ref config
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet_test.2024-02-26.tsv"

# Create output and log directories if they don't exist
input_fastq_dir="./fastq"
output_bam_dir="./bam"
log_dir="./bam/log"
mkdir -p ${output_bam_dir} && chmod a+rw ${output_bam_dir}
mkdir -p ${log_dir} && chmod a+rw ${log_dir}

############### function ###############
run_bwa(){
    Sample_ID=${1}
    
    bwa mem \
    -t ${threads} \
    -T 0 \
    -R "@RG\tID:${Sample_ID}\tSM:${Sample_ID}\tPL:ILLUMINA" \
    ${ref_fa} \
    "${input_fastq_dir}/${Sample_ID}_1.fq.gz" \
    "${input_fastq_dir}/${Sample_ID}_2.fq.gz" \
    | samtools view -Shb -o "${output_bam_dir}/${Sample_ID}.bam" - 
}

############### main ###############
# col 6 case_IDs TCGA-4G-AAZR sort 並去除重複
case_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# for sample sample_IDs
for case_IDs in "${case_IDs[@]}"
do
    echo "Processing ${case_IDs}"
    tumor="${case_IDs}-01A"
    normal="${case_IDs}-11A"
	
    run_bwa ${tumor} >> "${log_dir}/${tumor}.log" 2>&1 || { echo "Error processing ${tumor}"; exit 1; }
    run_bwa ${normal} >> "${log_dir}/${normal}.log" 2>&1 || { echo "Error processing ${normal}"; exit 1; }
done
