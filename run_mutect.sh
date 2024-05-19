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
germline="${db_dir}/af-only-gnomad.hg38.vcf.gz"
pon="${db_dir}/1000g_pon.hg38.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"
small_exac="${db_dir}/small_exac_common_3.hg38.vcf.gz"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet_test.2024-02-26.tsv"

# Create output and log directories if they don't exist
input_bam_dir="${project_dir}/sorted_du_bqsr_bam"
output_mutect_dir="./mutect"
log_dir="./mutect/log"
mkdir -p ${output_mutect_dir} && chmod a+rw ${output_mutect_dir}
mkdir -p ${log_dir}g && chmod a+rw ${log_dir}

############### function ###############
run_mutect(){
    case_ID=${1}
    tumor_file="${input_bam_dir}/${case_ID}-01A.sorted.du.bqsr.bam"
    normal_file="${input_bam_dir}/${case_ID}-11A.sorted.du.bqsr.bam"
    
    {
        gatk Mutect2 \
            -R ${ref_fa} \
            -I ${tumor_file} \
            -I ${normal_file} \
            -tumor "${case_ID}-01A" \
            -normal "${case_ID}-11A" \
            --germline-resource ${germline} \
            --native-pair-hmm-threads 20 \
            -pon ${pon} \
            -O "${output_mutect_dir}/${case_ID}.vcf"

        gatk GetPileupSummaries \
            -I ${tumor_file} \
            -V ${small_exac} \
            -L ${small_exac} \
            -O "${output_mutect_dir}/${case_ID}-01A.pileups.table"

        gatk GetPileupSummaries \
            -I ${normal_file} \
            -V ${small_exac} \
            -L ${small_exac} \
            -O "${output_mutect_dir}/${case_ID}-11A.pileups.table"

        gatk CalculateContamination \
            -I "${output_mutect_dir}/${case_ID}-01A.pileups.table" \
            -matched "${output_mutect_dir}/${case_ID}-11A.pileups.table" \
            -O "${output_mutect_dir}/${case_ID}.contamination.table" \
            -tumor-segmentation "${output_mutect_dir}/${case_ID}.segments.table"

        gatk FilterMutectCalls \
            -V "${output_mutect_dir}/${case_ID}.vcf" \
            --tumor-segmentation "${output_mutect_dir}/${case_ID}.segments.table" \
            --contamination-table "${output_mutect_dir}/${case_ID}.contamination.table" \
            -O "${output_mutect_dir}/${case_ID}.filtered.vcf" \
            -R ${ref_fa}

    } > "${log_dir}/${case_ID}.log" 2>&1 || {
        echo "Error processing ${case_ID}. Check log file for details."
    }
}

############### main ###############
# Extract sample_IDs and remove duplicates
sample_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# Run mutect for each sample
for sample_ID in "${case_IDs[@]}"
do  
    echo "Processing ${case_IDs}"
    run_mutect "${case_IDs}"
done
