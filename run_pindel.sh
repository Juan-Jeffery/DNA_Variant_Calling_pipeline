#!/bin/bash
# 20240529

############### preprocessing ###############
# config 
env_name="WES_Indel" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet.2024-02-26.tsv"

# Create output and log directories if they don't exist
pindel_dir="./variants_calling/pindel"
log_dir="${pindel_dir}/log"
mkdir -p ${pindel_dir}
mkdir -p ${log_dir}

############### function ###############
run_pindel(){
    case_ID=${1}
    tumor=${2}
    normal=${3}

    config_file="${pindel_dir}/${case_ID}.config.txt"
    echo "${bam_input_dir}/${normal}.bam 250 ${case_ID}_N" > ${config_file}
    echo "${bam_input_dir}/${tumor}.bam 250 ${case_ID}_T" >> ${config_file}
    
    {
        # step 1: pindel variant calling
        pindel -f ${ref_fa} -i ${config_file} -c ALL -T 28 -o ${pindel_dir}/${case_ID}.Pindel

        # step 2: extract indel summary lines
        grep ChrID ${pindel_dir}/${case_ID}.Pindel_D > ${pindel_dir}/${case_ID}.D.head
        grep ChrID ${pindel_dir}/${case_ID}.Pindel_SI > ${pindel_dir}/${case_ID}.SI.head
        cat ${pindel_dir}/${case_ID}.D.head ${pindel_dir}/${case_ID}.SI.head > ${pindel_dir}/${case_ID}.all.head

        # run pindel2vcf
        pindel2vcf -R ${ref_fa} -r ${ref_fa} -p ${pindel_dir}/${case_ID}.all.head -d 20240325 -v ${pindel_dir}/${case_ID}.d_SI.Pindel.vcf -G

    } > ${log_dir}/${case_ID}.log 2>&1 || { echo "Error processing ${case_ID}. Check log file for details."; exit 1; }
}

############### main ###############
# Extract case_IDs and remove duplicates
case_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# for sampleID in case_IDs
for case_ID in "${case_IDs[@]}"
do
    echo "Processing ${case_ID}..."
    tumor="${case_ID}-01A"
    normal="${case_ID}-11A"
    run_pindel ${case_ID} ${tumor} ${normal}
done

echo "All processes completed successfully."
