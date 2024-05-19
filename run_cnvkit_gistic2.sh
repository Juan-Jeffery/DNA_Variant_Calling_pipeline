#!/bin/bash
# 20240529

############### preprocessing ###############
# config 
env_name="WES_CNV" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

threads="20"

# ref config
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa" # 修改ref位置
gistic_dir="/path/to/gistic" # Ensure gistic_dir is correctly defined
refgene="${gistic_dir}/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"

# project config
project_dir="/home/data/dataset/CHOL_10sample" # 專案資料夾
cd ${project_dir}
access="./cnv_ref/access-10kb.hg38.bed"
refFlat="./cnv_ref/refFlat.txt"
illumina_bed="./cnv_ref/Twist_ILMN_Exome_2.0_Plus_Panel.hg38.bed"

# Create output and log directories if they don't exist
bam_input_dir="./sorted_du_bqsr_bam"
cnv_output_dir="./cnv"
gistic_output_dir="./gistic"
log_dir="./cnv/log"
mkdir -p ${cnv_output_dir}
mkdir -p ${gistic_output_dir}
mkdir -p ${log_dir}

############### function ###############
run_cnvkit(){
    {
        echo "Running CNVkit..."
        # call copy number variants
        cnvkit.py batch ${project_dir}/${bam_input_dir}/*01A*.bam \
            --normal ${project_dir}/${bam_input_dir}/*11A*.bam \
            --targets ${illumina_bed} --annotate ${refFlat} \
            --fasta ${ref_fa} --access ${access} -p ${threads} \
            --output-reference ${cnv_output_dir}/cnv_ref.cnn --output-dir ${cnv_output_dir} \
            --diagram --scatter
        
        echo "Exporting to GISTIC format..."
        # export to gistic format
        cnvkit.py export seg ${cnv_output_dir}/*.cns \
            -o ${cnv_output_dir}/gistic.segments
    } > ${log_dir}/cnvkit.log 2>&1 || { echo "Error in CNVkit"; exit 1; }
}

run_gistic2(){
    {
        echo "Running GISTIC2..."
        # run gistic
        gistic2 -b ${gistic_output_dir} \
            -seg ${cnv_output_dir}/gistic.segments \
            -refgene ${refgene} \
            -genegistic 1 \
            -smallmem 1 \
            -broad 1 \
            -brlen 0.5 \
            -conf 0.90 -armpeel 1 -savegene 1 -gcm "extreme"
    } > ${log_dir}/gistic2.log 2>&1 || { echo "Error in GISTIC2"; exit 1; }
}

############### main ###############
# call functions 
run_cnvkit
run_gistic2

echo "All processes completed successfully."
