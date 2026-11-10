library(data.table)
library(stringr)
library(arrow)
setwd("/data/cephfs-1/work/projects/sanders-digital-karyotype/pipelines/mosaicatcher-pipeline-BIH")
source("workflow/scripts/mosaiclassifier_scripts/mosaiClassifier/mosaiClassifier.R")


################################################################################
# inputs
################################################################################
args <- commandArgs(trailingOnly = TRUE)
rdata_file <- args[1]
par_out <- args[2]


read_save_parquet <- function(file, out_file) {
    dt <- readRDS(file)
    message(stringr::str_glue("RDS read..."))
    dt <- mosaiClassifierPostProcessing(dt, haplotypeMode=T, regularizationFactor=1e-10)
    
    # pop priors
    dt[, pop.prior := sum(nb_hap_ll), by = .(chrom, start, end, sample, haplotype)]
    dt[, pop.prior := pop.prior / sum(pop.prior), by = .(chrom, start, end, sample, cell)]
    dt[, nb_hap_pp := nb_hap_pp * pop.prior]
    print("pop_priors applied")

    # haplotag_probs
    dt[!is.na(haplotag.prob), nb_hap_pp := nb_hap_pp * haplotag.prob]
    print("haplotag probs applied")

    # normalising each prob dist
    dt[, nb_hap_pp := nb_hap_pp / sum(nb_hap_pp), by = .(sample, cell, chrom, start, end)]

    dt <- dt[, .(chrom, start, end, sample, cell, num_bins, haplotype, haplotag.prob, nb_hap_pp)]

    message(stringr::str_glue("Saving to file: {out_file}"))
    arrow::write_parquet(x = dt,
                         sink = out_file,
                         compression = "brotli",
                         use_dictionary = FALSE
                         )
}

read_save_parquet(rdata_file, par_out)
