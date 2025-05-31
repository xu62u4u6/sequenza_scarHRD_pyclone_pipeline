#%%
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
# nextflow dir
nextflow_dir = Path.cwd().parent
project_dir = "/home/data/data_dingyangliu/ASC_0217/"
# load data
os.chdir(project_dir)
from pymaftools.pymaftools import * 
from scripts.utils.setup_env import ASCDataLoader
dataloader = ASCDataLoader(project_dir)

all_sample_metadata = dataloader.all_sample_metadata
AST_sample_metadata = dataloader.AST_sample_metadata
cnv_table = dataloader.cnv_table
snv_table = dataloader.snv_table

snv_ATS_subset, cnv_ATS_subset = dataloader.get_ATS_subset()
#%%
def infer_normal_cn(row):
    chrom = row["chromosome"]
    sex = row["sex"]

    if chrom in [f"chr{i}" for i in range(1, 23)]:
        return 2
    elif chrom == "chrX":
        if sex == "M":
            return 1
        elif sex == "F":
            return 2
    elif chrom == "chrY":
        if sex == "M":
            return 1
        elif sex == "F":    
            return 0
    else:
        return np.nan


def prepare_pyclone_input(sequenza_df: pd.DataFrame, snp_df: pd.DataFrame) -> pd.DataFrame:
    """
    給定 Sequenza CNV 結果與 SNP 資料，整合成符合 PyClone-VI 格式的輸入。
    
    Parameters:
    - sequenza_df: Sequenza 的 segment 資料（包含 CNt, A, B, normal_cn）
    - snp_df: MAF 或其他格式的 SNP 資料（包含 ref/alt count）

    Returns:
    - DataFrame: 包含 mutation_id、ref_counts、var_counts、normal_cn、major/minor_cn、sample_ID 的 PyClone-VI 格式
    """

    sequenza_cnv = sequenza_merged_df[[
        "chromosome", "start.pos", "end.pos", "sample.id", "CNt", "A", "B", "normal_cn", "cellularity"
    ]].rename(columns={
        "chromosome": "Chromosome",
        "start.pos": "start",
        "end.pos": "end",
        "sample.id": "sample_ID",
        "CNt": "total_cn",
        "A": "major_cn",
        "B": "minor_cn",
        "cellularity": "tumour_content"
    })

    snp_slim = snp_df.reset_index()[[
        "Hugo_Symbol", "Chromosome", "Start_Position",
        "Reference_Allele", "Tumor_Seq_Allele2",
        "t_ref_count", "t_alt_count", "sample_ID"
    ]].copy()

    # Step 2: 染色體排序定義
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chrom_dtype = pd.api.types.CategoricalDtype(categories=chrom_order, ordered=True)

    # 強制轉換排序類別
    snp_slim["Chromosome"] = snp_slim["Chromosome"].astype(chrom_dtype)
    sequenza_cnv["Chromosome"] = sequenza_cnv["Chromosome"].astype(chrom_dtype)

    # Step 3: 合併資料（SNP 對應 CNV）
    results = []
    for sample_id in snp_slim["sample_ID"].unique():
        snp_sample = snp_slim[snp_slim["sample_ID"] == sample_id]
        cnv_sample = sequenza_cnv[sequenza_cnv["sample_ID"] == sample_id]

        for chrom in snp_sample["Chromosome"].dropna().unique():
            snp_chrom = snp_sample[snp_sample["Chromosome"] == chrom].copy()
            cnv_chrom = cnv_sample[cnv_sample["Chromosome"] == chrom].copy()

            if snp_chrom.empty or cnv_chrom.empty:
                continue

            snp_chrom = snp_chrom.sort_values("Start_Position")
            cnv_chrom = cnv_chrom.sort_values("start")

            matched = pd.merge_asof(
                snp_chrom,
                cnv_chrom,
                left_on="Start_Position",
                right_on="start",
                direction="backward",
                allow_exact_matches=True,
                suffixes=('', '_cnv')
            )

            # 過濾不落在區段內的 SNP
            matched = matched[matched["Start_Position"] <= matched["end"]]
            results.append(matched)


    matched_df = pd.concat(results, ignore_index=True)

    # Step 4: 建立 mutation_id 與輸出格式
    matched_df["mutation_id"] = (
        matched_df["Hugo_Symbol"] + ":" +
        matched_df["Chromosome"].astype(str) + ":" +
        matched_df["Start_Position"].astype(str) + ":" +
        matched_df["Reference_Allele"] + ">" +
        matched_df["Tumor_Seq_Allele2"]
    )

    pyclone_input = matched_df[[
        "mutation_id",
        "sample_ID",       # 要改成 sample_id
        "t_ref_count",     # 要改成 ref_counts
        "t_alt_count",     # 要改成 alt_counts
        "major_cn", "minor_cn", "normal_cn",
        "tumour_content"      # optional
    ]].rename(columns={
        "t_ref_count": "ref_counts",
        "t_alt_count": "alt_counts",
        "sample_ID": "sample_id",
    })

    return pyclone_input
#%%
# read MAF
maf = MAF.read_csv("data/WES/all_case_maf.maf")
snp = maf[maf.Variant_Type == "SNP"]
all_sample_metadata = all_sample_metadata.set_index("case_ID")
AST_case_ID = all_sample_metadata.loc[all_sample_metadata.AS_sample == True].index

# preprocessing snv data
snp["case_ID"] = snp["sample_ID"].str.rsplit("_", n=1).str[0]
snp = snp.loc[snp.case_ID.isin(all_sample_metadata.index)] # drop subtype B

#%%
sequenza_path = os.path.join(nextflow_dir, "output", "all_segments_merged.tsv")
sequenza_merged_df = pd.read_csv(sequenza_path, sep="\t")
sequenza_merged_df["case_ID"] = sequenza_merged_df["sample.id"].str.rsplit("_", n=1).str[0]
sequenza_merged_df = sequenza_merged_df.loc[sequenza_merged_df.case_ID.isin(all_sample_metadata.index)]
sequenza_merged_df["sex"] = sequenza_merged_df["case_ID"].apply(lambda case_ID: all_sample_metadata.loc[case_ID, "sex"])
sequenza_merged_df["normal_cn"] = sequenza_merged_df.apply(infer_normal_cn, axis=1)

#%%
pyclone_input_df = prepare_pyclone_input(sequenza_merged_df, snp)
pyclone_input_df = pyclone_input_df.loc[pyclone_input_df.major_cn != 0]
for col in ["major_cn", "minor_cn", "normal_cn"]:
    pyclone_input_df[col] = pyclone_input_df[col].astype(int)
# pyclone_input_df.to_csv("test_pyclone/pyclone_input.tsv", sep="\t", index=False)

pyclone_input_df["case_ID"] = pyclone_input_df["sample_id"].str.rsplit("_", n=1).str[0]

# make dir
os.makedirs("data", exist_ok=True)
os.makedirs("data/pyclone_input", exist_ok=True)

for case_ID in AST_case_ID:
    sub_df = pyclone_input_df.loc[pyclone_input_df["case_ID"] == case_ID]
    sub_df.to_csv(f"data/pyclone_input/{case_ID}.tsv", sep="\t", index=False)