#!/bin/bash
#20240529
############### preprocessing ############### 
# config 
env_name="WES_preprocessing_SNP" # 環境名稱請確認
# enable bash script calling conda function
eval "$(conda shell.bash hook)"
conda activate ${env_name}

threads=20
picard=$(find "/home/${USER}/miniconda3/envs/${env_name}/share/" -name "picard.jar" | head -n 1) # 找到picard.jar

# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"
ref_gtf="/home/data/dataset/BWA_index/gencode.v36.annotation.gtf"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet_test.2024-02-26.tsv"

# Create output and log directories if they don't exist
input_bam_dir="./bam"
output_bam_dir="./sorted_du_bqsr_bam"
log_dir="./sorted_du_bqsr_bam/log"
mkdir -p ${output_bam_dir} && chmod a+rw ${output_bam_dir}
mkdir -p ${log_dir} && chmod a+rw ${log_dir}

############### function ###############
bqsr(){
    # input: .bam
    # intermediate: .sorted.bam, .sorted.du.bam, .sorted.du.matrix.txt, recalibration_report.txt
    # output: .sorted.du.bqsr.bam

    sample_ID=${1}

    # step1. bam file排序
    java -Xmx50g -jar ${picard} SortSam \
        I="${input_bam_dir}/${sample_ID}.bam" \
        O="${output_bam_dir}/${sample_ID}.sorted.bam" \
        SO=coordinate \
        VALIDATION_STRINGENCY=STRICT \
        CREATE_INDEX=true
        
    # mapping rate (QC)
    samtools flagstat "${output_bam_dir}/${sample_ID}.sorted.bam" "${output_bam_dir}/${sample_ID}.mapping_rate_human.txt"

    # step2. 將重複序列進行標注
    java -Xmx50g -jar ${picard} MarkDuplicates \
        I="${output_bam_dir}/${sample_ID}.sorted.bam" \
        O="${output_bam_dir}/${sample_ID}.sorted.du.bam" \
        M="${output_bam_dir}/${sample_ID}.sorted.du.matrix.txt" \
        REMOVE_DUPLICATES=false \
        CREATE_INDEX=true 2>> "${log_dir}/MARK_DUP.log"

    # step3. 校正base品質
    gatk BaseRecalibrator \
        -R ${ref_fa} \
        -I "${output_bam_dir}/${sample_ID}.sorted.du.bam" \
        -O "${output_bam_dir}/${sample_ID}.recalibration_report.txt" \
        --known-sites ${dbsnp}

    # step4. 套用校正
    gatk ApplyBQSR \
        -R ${ref_fa} \
        -I "${output_bam_dir}/${sample_ID}.sorted.du.bam" \
        -O "${output_bam_dir}/${sample_ID}.sorted.du.bqsr.bam" \
        --bqsr-recal-file "${output_bam_dir}/${sample_ID}.recalibration_report.txt"
    
    # Reads QC
    qualimap bamqc \
        -bam "${output_bam_dir}/${sample_ID}.sorted.du.bqsr.bam" \ 
        -gff ${ref_gtf} \
        -outfile "${output_bam_dir}/${sample_ID}.qualimap_sample.pdf" \
        -sd \
        -nt 20 \
        --java-mem-size=64G \
}

############### main ###############
# col 6 sample_ID TCGA-4G-AAZR sort 並去除重複
sample_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# for each sample_ID
for sample_ID in "${sample_IDs[@]}"
do
    echo "Processing ${sample_ID}"
    tumor="${sample_ID}-01A"
    normal="${sample_ID}-11A"
    
    bqsr ${tumor} >> "${log_dir}/${tumor}.log" 2>&1 || { echo "Error processing ${tumor}"; exit 1; }
    bqsr ${normal} >> "${log_dir}/${normal}.log" 2>&1 || { echo "Error processing ${normal}"; exit 1; }
done
