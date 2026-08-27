library(data.table)
library(stringr)
library(arrow)
# setwd("/data/cephfs-1/work/projects/sanders-digital-karyotype/mosaicatcher-pipeline-BIH")
source("workflow/scripts/mosaiclassifier_scripts/mosaiClassifier/mosaiClassifier.R")


# data.table::setDTthreads(15)
# data.table::getDTthreads()
################################################################################
# inputs
################################################################################
args <- commandArgs(trailingOnly = TRUE)
rdata_file <- args[1]
par_out <- args[2]
# rdata_file <- "~/sanders-digital-karyotype/Digital-Karyotype-SubCohorts/RPE-TALL"


read_save_parquet <- function(file, out_file) {
    dt <- readRDS(file)
    message(stringr::str_glue("RDS read..."))
    dt <- mosaiClassifierPostProcessing(dt, 1e-10)

    # keeping only useful columns
    dt <- dt[, .(chrom, start, end, sample, cell, num_bins, haplotype, nb_hap_pp)]

    sample <- dt[, unique(sample)]
    # out_file <- file.path(dirname(file),
    #                       "probabilities.parquet"
    # )
    message(stringr::str_glue("Saving to file: {out_file}"))
    # data.table::fwrite(dt, out_file, sep = "\t")
    arrow::write_parquet(x = dt,
                         sink = out_file,
                         compression = "brotli",
                         use_dictionary = FALSE
    )
}

read_save_parquet(rdata_file, par_out)
