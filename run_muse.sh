#!/bin/bash
#20240529
############### preprocessing ###############
# config 
env_name="WES_preprocessing_SNP" # 環境名稱請確認
# enable bash script calling conda function
eval "$(conda shell.bash hook)"
conda activate ${env_name}

# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet_test.2024-02-26.tsv"


# Create output and log directories if they don't exist
input_bam_dir="${project_dir}/sorted_du_bqsr_bam"
output_muse_dir="./MuSE"
log_dir="./MuSE/log"
mkdir -p ${muse_dir} && chmod a+rw ${muse_dir}
mkdir -p ${log_dir} && chmod a+rw ${log_dir}
############### function ###############
run_muse(){
    sample_ID=${1}

    tumor_file="${bam_dir}/${case_IDs}-01A.sorted.du.bqsr.bam"
    normal_file="${bam_dir}/${case_IDs}-11A.sorted.du.bqsr.bam"

    MuSE call -f ${ref_fa} \
    -O ${muse_dir}/${sample_ID} \
    ${tumor_file} ${normal_file} \
    > "${muse_dir}/log/${case_IDs}.call.log" 2>&1

    MuSE_exit_status=$?
    if [ "${MuSE_exit_status}" -eq 0 ]; then
        echo "${case_IDs} MuSE call completed successfully" >> "${muse_dir}/muse.log"
    else
        echo "${case_IDs} MuSE call failed with exit status ${MuSE_exit_status}" >> "${muse_dir}/muse.log"
        return ${MuSE_exit_status}
    fi

    MuSE sump -I ${muse_dir}/${sample_ID}.MuSE.txt -E \
    -D ${dbsnp} \
    -O ${muse_dir}/${sample_ID}.MuSE.vcf \
    >> "${muse_dir}/log/${case_IDs}.sump.log" 2>&1

    MuSE_exit_status=$?
    if [ "${MuSE_exit_status}" -eq 0 ]; then
        echo "${case_IDs} MuSE sump completed successfully" >> "${muse_dir}/muse.log"
    else
        echo "${case_IDs} MuSE sump failed with exit status ${MuSE_exit_status}" >> "${muse_dir}/muse.log"
    fi
}

############### main ###############
# col 6 sample_ID TCGA-4G-AAZR sort並去除重複
sample_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# cp files
for case_IDs in "${case_IDs[@]}"
do
    normal="${case_IDs}-01A"
    tumor="${case_IDs}-11A"
    cp -n "${bam_dir}/${normal}.sorted.du.bqsr.bai" "${bam_dir}/${normal}.sorted.du.bqsr.bam.bai" 2>> "${muse_dir}/log/${case_IDs}.cp.log"
    cp -n "${bam_dir}/${tumor}.sorted.du.bqsr.bai" "${bam_dir}/${tumor}.sorted.du.bqsr.bam.bai" 2>> "${muse_dir}/log/${case_IDs}.cp.log"
done

# run muse
for case_IDs in "${case_IDs[@]}"
do  
    echo "Starting MuSE for ${case_IDs}" >> "${muse_dir}/muse.log"
    run_muse "${case_IDs}"
done

# Wait for all background processes to finish
wait
echo "All MuSE jobs completed." >> "${muse_dir}/muse.log"