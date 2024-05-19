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
output_QC_dir="./QC"
mkdir -p ${output_QC_dir} && chmod a+rw ${output_QC_dir}

############### function ###############
run_fastqc(){
    Sample_id=${tumor}
    fastqc -o ${output_QC_dir} ${input_fastq_dir}/${Sample_id}.fastq.gz
}

run_multiqc(){
    cd ${output_QC_dir}
    multiqc .
}
############### main ###############
# col 6 case_IDs TCGA-4G-AAZR sort 並去除重複(cut是從1開始算)
case_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# for sample sample_IDs
for case_IDs in "${case_IDs[@]}"
do
    echo "Processing ${case_IDs}"
    tumor="${case_IDs}-01A"
    normal="${case_IDs}-11A"
	
    run_fastqc ${tumor}
    run_fastqc ${normal}
done

run_multiqc


