#!/bin/bash
#20240529
############### preprocessing ###############
# config 
env_name="WES_preprocessing_SNP" # 環境名稱請確認
# enable bash script calling conda function
eval "$(conda shell.bash hook)"
conda activate ${env_name}

# ref
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa" # 修改ref位置
varscan="/home/${USER}/miniconda3/envs/${env_name}/share/varscan-*/VarScan.jar" # 抓取使用者的環境資料夾

# project config
project_dir="/home/data/dataset/CHOL_10sample" # 專案資料夾
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet_test.2024-02-26.tsv"

# Create output and log directories if they don't exist
input_bam_dir="./sorted_du_bqsr_bam_biobambam"
output_varscan_dir="./variants_calling/varscan"
log_dir="./variants_calling/varscan/log"
mkdir -p ${output_varscan_dir} && chmod a+rw ${output_varscan_dir}
mkdir -p ${log_dir} && chmod a+rw ${log_dir}

############### function ###############
run_varscan(){
    case_ID=${1}
    tumor="${case_ID}-01A"
    normal="${case_ID}-11A"

    # Step 1: Mpileup; Samtools
    samtools mpileup \
    -f ${ref_fa} \
    -q 1 \
    -B  \
    ${input_bam_dir}/${normal}.sorted.du.bqsr.bam \
    ${input_bam_dir}/${tumor}.sorted.du.bqsr.bam \
    > ${output_varscan_dir}/${case_ID}_intermediate_mpileup.pileup

    echo "step1 samtools mpileup ${case_ID} done" >> ${log_dir}/${case_ID}.log
    
    # Step 2: Varscan Somatic; Varscan.v2
    java -jar ${varscan} somatic \
    "${output_varscan_dir}/${case_ID}_intermediate_mpileup.pileup" \
    "${output_varscan_dir}/${case_ID}.varscan" \
    --mpileup 1 \
    --min-var-freq 0.01 \
    --output-vcf >> ${log_dir}/${case_ID}.log 2>&1
    echo "step2 somatic ${case_ID} done" >> ${log_dir}/${case_ID}.log
    
    # Step 3: Varscan ProcessSomatic; Varscan.v2
    java -jar ${varscan} processSomatic \
    "${output_varscan_dir}/${case_ID}.varscan.snp.vcf" \
    --min-tumor-freq 0.01 >> ${log_dir}/${case_ID}.log 2>&1
    
    java -jar ${varscan} processSomatic \
    "${output_varscan_dir}/${case_ID}.varscan.indel.vcf" \
    --min-tumor-freq 0.01 >> ${log_dir}/${case_ID}.log 2>&1
    echo "step3 processSomatic ${case_ID} done" >> ${log_dir}/${case_ID}.log

    # Step 4: del tmp
    rm "${output_varscan_dir}/${case_ID}_intermediate_mpileup.pileup"
    echo "step4 clean ${case_ID} tmp file done" >> ${log_dir}/${case_ID}.log

    # Step 5: Optional Concatenation of VCF files
    # Uncomment if needed
    # varscan_indel="${output_varscan_dir}/${case_ID}.varscan.indel.Somatic.hc.vcf"
    # varscan_snp="${output_varscan_dir}/${case_ID}.varscan.snp.Somatic.hc.vcf"
    # bcftools concat -o ${output_varscan_dir}/${case_ID}.varscan.vcf ${varscan_indel} ${varscan_snp}
}

############### main ###############
# load sample sheet 
case_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# for sampleID in case_IDs
for case_ID in "${case_IDs[@]}"
do
    run_varscan ${case_ID}
done

# Wait for all background jobs to finish
wait
echo "All VarScan jobs completed."
