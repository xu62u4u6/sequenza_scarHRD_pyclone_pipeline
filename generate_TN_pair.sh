#!/bin/bash
bam_dir="/home/data/data_dingyangliu/ASC_0217/bam"

# 輸出表頭
echo -e "sample_id\ttumor_bam\tnormal_bam"

# 遍歷所有子資料夾
find "$bam_dir" -mindepth 1 -maxdepth 1 -type d | while read case_dir; do
    case_id=$(basename "$case_dir")
    normal_bam="$case_dir/${case_id}_N.bam"
    # 若 normal 不存在，跳過
    [[ ! -f "$normal_bam" ]] && continue

    # 找出非 N 的 tumor bams
    for tumor_bam in "$case_dir"/${case_id}_[AST].bam; do
        [[ ! -e "$tumor_bam" ]] && continue

        tumor_id=$(basename "$tumor_bam" .bam)
        echo -e "${tumor_id}\t${tumor_bam}\t${normal_bam}"
    done
done
