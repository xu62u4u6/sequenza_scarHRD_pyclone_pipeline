# config
options(future.globals.maxSize= 100e+9)
Sys.setenv("VROOM_CONNECTION_SIZE"=999999999)  # GB 1,000,000,000
args <- commandArgs(trailingOnly=TRUE)
sample_id <- args[1]
seqz <- args[2]
output_dir <- "scarHRD_result"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ====== Load package ======
library("scarHRD")

# ====== scarHRD ======
scar_score(seqz, reference = "grch38", seqz=TRUE, outputdir=output_dir)