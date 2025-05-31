# config
options(future.globals.maxSize = 100e+9)
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999999)
args <- commandArgs(trailingOnly=TRUE)
sample_id <- args[1]
input_file <- args[2]
n_cores <- as.integer(args[3])
output_dir <- "sequenza_result"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
library(sequenza)
library(data.table)

# Step 1: Extract
seqz_data <- sequenza.extract(input_file, chromosome.list = paste0("chr", c(1:22, "X", "Y")), parallel = n_cores)

# Step 2: Fit
cp_table <- sequenza.fit(seqz_data, mc.cores = n_cores)
cint <- get.ci(cp_table)
cellularity <- cint$max.cellularity
ploidy <- cint$max.ploidy

# Step 3: Output
sequenza.results(
    sequenza.extract = seqz_data,
    cp.table = cp_table,
    sample.id = sample_id,
    out.dir = output_dir
)

# Step 4: Modify segments
segments_file <- file.path(output_dir, paste0(sample_id, "_segments.txt"))
if (!file.exists(segments_file)) {
    stop(paste("Segments file does not exist:", segments_file))
}
segments <- fread(segments_file)
segments$cellularity <- cellularity
segments$ploidy <- ploidy
segments$sample.id <- sample_id

output_file <- file.path(output_dir, paste0(sample_id, "_segments_modified.txt"))
fwrite(segments, output_file, sep = "\t")