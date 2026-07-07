#!/usr/bin/Rscript
options(error = traceback)
args <- commandArgs(TRUE)

# ensure lib path / packages resolved
.libPaths(c(.libPaths(), args[6]))

# load package (will stop if unavailable)
suppressPackageStartupMessages({
  if (!requireNamespace("StrandPhaseR", quietly = TRUE)) {
    stop("StrandPhaseR package is not available in the R library paths. Check .libPaths and args[6].")
  }
  library(StrandPhaseR)
})

# sanity log for Snakemake
message("Session info (first lines):")
message(paste(head(capture.output(sessionInfo()), 20), collapse = "\n"))
message("Loaded packages: ", paste(.packages(), collapse = ", "))

# define the patched function locally (exact same body as before)
patched_exportConsensus <- function (data.bases, data.quals, min.cov = 2, translateBases = FALSE) 
{
  indices <- which(data.bases != 0, arr.ind = TRUE)
  values <- data.bases[indices]
  col.vals <- split(values, (indices[, 2]))
  
  indices <- which(data.quals != 0, arr.ind = TRUE)
  values <- data.quals[indices]
  col.quals <- split(values, (indices[, 2]))
  
  base.freq <- lapply(col.vals, table)
  positions <- as.numeric(names(base.freq))
  max.cov <- sapply(base.freq, function(x) max(x))
  score <- sapply(base.freq, function(x) sum(x) - max(x))
  
  mask <- which(max.cov > score & max.cov >= min.cov)
  
  base.freq <- base.freq[names(mask)]
  col.quals <- col.quals[names(mask)]
  col.vals <- col.vals[names(mask)]
  positions <- as.numeric(names(mask))
  
  if (length(positions) == 0) {
    return(0)
  }
  
  entropy <- sapply(col.vals, calcEnt)
  
  # ---- patched & logged scoring block ----
  out_list <- mapply(calcProb, bases = col.vals, quals = col.quals, SIMPLIFY = FALSE)
  
  keep <- which(vapply(out_list, function(x) length(x[[2]]), integer(1)) == 1L)
  
  if (length(keep) < length(out_list)) {
    dropped <- length(out_list) - length(keep)
    dropped_pos <- as.numeric(names(out_list)[setdiff(seq_along(out_list), keep)])
    
    message(
      sprintf(
        "exportConsensus: dropping %d unscorable sites (example positions: %s)",
        dropped,
        paste(head(dropped_pos, 3), collapse = ",")
      )
    )
  }
  
  # subset consistently
  positions <- positions[keep]
  entropy   <- entropy[keep]
  col.vals  <- col.vals[keep]
  
  cov <- sapply(col.vals, length)
  
  bases  <- vapply(out_list[keep], function(x) x[[1]], integer(1))
  scores <- vapply(out_list[keep], function(x) x[[2]], numeric(1))
  
  if (translateBases) {
    bases <- chartr("1234", "ACGT", as.character(bases))
  }
  
  assem.haps <- data.frame(
    pos   = positions,
    bases = bases,
    cov   = cov,
    score = scores,
    ent   = entropy
  )
  
  rownames(assem.haps) <- NULL
  assem.haps
}

# try to insert the patched function into the package namespace
ns_name <- "StrandPhaseR"
ns <- asNamespace(ns_name)

ok <- FALSE
try({
  unlockBinding("exportConsensus", ns)
  assign("exportConsensus", patched_exportConsensus, envir = ns)
  lockBinding("exportConsensus", ns)
  message("Patched exportConsensus successfully in namespace '", ns_name, "'.")
}, silent = FALSE)

try({
  ok <- identical(get("exportConsensus", envir = ns), patched_exportConsensus)
}, silent = TRUE)

message("Verification that exportConsensus was replaced: ", ok)
stopifnot(ok)  # enable temporarily to ensure it's really patched

# -------------------------
# Patched phaseChromosome()
# -------------------------
patched_phaseChromosome <- function (inputfolder, outputfolder = "./StrandPhaseR_analysis", 
    positions = NULL, WCregions = NULL, chromosome = NULL, pairedEndReads = TRUE, 
    min.mapq = 10, min.baseq = 20, num.iterations = 2, translateBases = TRUE, 
    concordance = 0.9, fillMissAllele = NULL, splitPhasedReads = FALSE, 
    compareSingleCells = FALSE, exportVCF = NULL, bsGenome = NULL, 
    ref.fasta = NULL, assume.biallelic = FALSE) 
{
    message("Working on chromosome ", chromosome)

    ## --- BEGIN: seqlevels / seqlengths hardening patch ---
    # Reduce GRanges to the single chromosome (safe pruning)
    if (!is.null(positions)) {
      positions <- keepSeqlevels(positions, chromosome, pruning.mode = "coarse")
    }
    if (!is.null(WCregions)) {
      WCregions <- keepSeqlevels(WCregions, chromosome, pruning.mode = "coarse")
    }

    # If seqlength for this chromosome is missing (NA), try to extract from first BAM
    pos_sl <- if (!is.null(positions)) seqlengths(positions)[as.character(chromosome)] else NA
    wc_sl  <- if (!is.null(WCregions)) seqlengths(WCregions)[as.character(chromosome)] else NA

    if (is.na(pos_sl) || is.na(wc_sl)) {
      # attempt to read BAM header referenced in WCregions to get chr length
      if (!is.null(WCregions) && length(WCregions) > 0 && !is.null(mcols(WCregions)$filename)) {
        first_bam_name <- as.character(mcols(WCregions)$filename[1])
        bam_candidate1 <- first_bam_name
        bam_candidate2 <- file.path(inputfolder, first_bam_name)
        bam_path <- if (file.exists(bam_candidate1)) bam_candidate1 else if (file.exists(bam_candidate2)) bam_candidate2 else NA

        if (!is.na(bam_path)) {
          hdr_targets <- tryCatch(Rsamtools::scanBamHeader(bam_path)[[1]]$targets, error = function(e) NULL)
          if (!is.null(hdr_targets) && !is.na(hdr_targets[chromosome])) {
            chr_len <- as.integer(hdr_targets[chromosome])
            sl_vec <- setNames(chr_len, chromosome)
            if (!is.null(positions)) seqlengths(positions) <- sl_vec
            if (!is.null(WCregions)) seqlengths(WCregions) <- sl_vec
            message("Assigned seqlengths for ", chromosome, " from BAM header (", bam_path, "): ", chr_len)
          } else {
            warning("Unable to determine seqlength for ", chromosome, " from BAM header; continuing with NA seqlengths.")
          }
        } else {
          warning("Referenced BAM not found (tried '", bam_candidate1, "' and '", bam_candidate2, "'). If seqlengths are needed, provide a bsGenome or BAM accessible to the pipeline.")
        }
      } else {
        warning("WCregions has no filename metadata or is empty; cannot derive seqlengths automatically.")
      }
    }
    ## --- END: seqlevels / seqlengths hardening patch ---

    phased.store <- file.path(outputfolder, "Phased")
    if (!dir.exists(phased.store)) {
        dir.create(phased.store, recursive = TRUE)
    }
    data.store <- file.path(outputfolder, "data")
    if (!dir.exists(data.store)) {
        dir.create(data.store, recursive = TRUE)
    }
    browser.store <- file.path(outputfolder, "browserFiles")
    if (!dir.exists(browser.store)) {
        dir.create(browser.store, recursive = TRUE)
    }
    vcf.store <- file.path(outputfolder, "VCFfiles")
    if (!dir.exists(vcf.store)) {
        dir.create(vcf.store, recursive = TRUE)
    }
    singlecell.store <- file.path(outputfolder, "SingleCellHaps")
    if (!dir.exists(singlecell.store)) {
        dir.create(singlecell.store, recursive = TRUE)
    }

    matrices <- loadMatrices(inputfolder = inputfolder, positions = positions, 
        WCregions = WCregions, pairedEndReads = pairedEndReads, 
        min.mapq = min.mapq, min.baseq = min.baseq)
    if (length(matrices) > 0) {
        srt.matrices <- sortMatrices(data.object = matrices, 
            num.iterations = num.iterations)
        assem.haps <- assembleHaps(data.object = srt.matrices, 
            translateBases = translateBases, concordance = concordance)
        if (!is.null(fillMissAllele)) {
            header <- utils::read.table(fillMissAllele, stringsAsFactors = FALSE, 
                fill = TRUE, comment.char = "&", nrows = 1)
            if (grepl(header, pattern = "VCF", ignore.case = TRUE)) {
                assem.haps <- fillGapsWithVCF(data.object = assem.haps, 
                  ref.vcf = fillMissAllele, chromosome = chromosome)
            }
            if (grepl(fillMissAllele, pattern = "\\.bam$")) {
                assem.haps <- fillGapsWithBam(data.object = assem.haps, 
                  merged.bam = fillMissAllele, min.mapq = min.mapq, 
                  min.baseq = min.baseq, translateBases = translateBases, 
                  chromosome = chromosome)
            }
        }
        if (compareSingleCells) {
            suppressWarnings(cell.comparisons.l <- compareSingleCellHaps(consensusHaps = assem.haps, 
                sortedHaps = srt.matrices, bin.size = 5))
            if (!is.null(cell.comparisons.l)) {
                destination <- file.path(singlecell.store, paste0(chromosome, 
                  "_singleCellHaps.pdf"))
                suppressWarnings(plotSingleCellHaps(data = cell.comparisons.l, 
                  file = destination))
                LOH.regions <- LOHseeker(data.object = cell.comparisons.l, 
                  chromosome = chromosome, bin.size = 5)
                LOH.regions.df <- data.frame(LOH.regions)
                destination <- file.path(singlecell.store, paste0(chromosome, 
                  "_singleCell_LOH.txt"))
                write.table(LOH.regions.df, file = destination, 
                  quote = F, row.names = F)
            }
        }
        chrName.hap1 <- data.frame(chr = rep(chromosome, nrow(assem.haps$hap1.cons)))
        assem.haps$hap1.cons <- cbind(chrName.hap1, assem.haps$hap1.cons)
        chrName.hap2 <- data.frame(chr = rep(chromosome, nrow(assem.haps$hap2.cons)))
        assem.haps$hap2.cons <- cbind(chrName.hap2, assem.haps$hap2.cons)
        destination <- file.path(data.store, paste0(chromosome, 
            "_phased.RData"))
        save(srt.matrices, file = destination)
        destination <- file.path(phased.store, paste0(chromosome, 
            "_phased_hap1.txt"))
        utils::write.table(assem.haps$hap1.cons, file = destination, 
            row.names = F)
        destination <- file.path(phased.store, paste0(chromosome, 
            "_phased_hap2.txt"))
        utils::write.table(assem.haps$hap2.cons, file = destination, 
            row.names = F)
        destination <- file.path(phased.store, paste0(chromosome, 
            "_phasedFiles_hap1.txt"))
        hap1.files <- data.frame(names(assem.haps$hap1.files), 
            do.call(rbind, lapply(assem.haps$hap1.files, rbind)))
        names(hap1.files) <- c("Filenames", "Simil", "Disimil")
        utils::write.table(hap1.files, file = destination, row.names = F)
        destination <- file.path(phased.store, paste0(chromosome, 
            "_phasedFiles_hap2.txt"))
        hap2.files <- data.frame(names(assem.haps$hap2.files), 
            do.call(rbind, lapply(assem.haps$hap2.files, rbind)))
        names(hap2.files) <- c("Filenames", "Simil", "Disimil")
        utils::write.table(hap2.files, file = destination, row.names = F)
        destination <- file.path(phased.store, "phased_haps.txt")
        utils::write.table(data.frame(assem.haps$assem.haps), 
            file = destination, row.names = F, col.names = F, 
            quote = F, append = T, sep = "\t")
        if (!is.null(exportVCF) & !is.null(bsGenome)) {
            exportVCF(index = exportVCF, outputfolder = vcf.store, 
                phasedHap = assem.haps, bsGenome = bsGenome, 
                chromosome = chromosome, assume.biallelic = assume.biallelic)
        }
        else if (!is.null(exportVCF) & !is.null(ref.fasta) & 
            is.null(bsGenome)) {
            exportVCF(index = exportVCF, outputfolder = vcf.store, 
                phasedHap = assem.haps, ref.fasta = ref.fasta, 
                chromosome = chromosome, assume.biallelic = assume.biallelic)
        }
        else if (!is.null(exportVCF) & is.null(bsGenome)) {
            exportVCF(index = exportVCF, outputfolder = vcf.store, 
                phasedHap = assem.haps, positions = positions, 
                chromosome = chromosome, assume.biallelic = assume.biallelic)
        }
        if (splitPhasedReads) {
            haps.gr <- splitReads(data.object = assem.haps, inputfolder = inputfolder, 
                pairedEndReads = pairedEndReads, min.mapq = 10, 
                filterAltAlign = TRUE)
            destination <- file.path(data.store, paste0(chromosome, 
                "_reads.RData"))
            save(haps.gr, file = destination)
            exportBedGraph(index = paste0(chromosome, "_hap1"), 
                outputfolder = browser.store, fragments = haps.gr$hap1, 
                col = "0,128,255")
            exportBedGraph(index = paste0(chromosome, "_hap2"), 
                outputfolder = browser.store, fragments = haps.gr$hap2, 
                col = "0,255,255")
        }
    }
    else {
        message(" Insufficient data to assemble haplotypes, skipping ...")
        if (!is.null(exportVCF)) {
            message("    Printing empty VCF file !!!")
            exportVCF(index = exportVCF, outputfolder = vcf.store, 
                positions = positions, bsGenome = bsGenome, chromosome = chromosome)
        }
    }
}

# -------------------------
# Inject into package namespace
# -------------------------
ns_name <- "StrandPhaseR"
ns <- asNamespace(ns_name)
ok <- FALSE

try({
  if (bindingIsLocked("phaseChromosome", ns)) {
    unlockBinding("phaseChromosome", ns)
  }
  assign("phaseChromosome", patched_phaseChromosome, envir = ns)
  lockBinding("phaseChromosome", ns)
  message("Patched phaseChromosome successfully in namespace '", ns_name, "'.")
}, silent = FALSE)

try({
  ok <- identical(get("phaseChromosome", envir = ns), patched_phaseChromosome)
}, silent = TRUE)

message("Verification that phaseChromosome was replaced: ", ok)
stopifnot(ok)  # fail fast if not patched

strandPhaseR(inputfolder = args[1], outputfolder = args[2], configfile = args[3], WCregions = args[4], positions = args[5], fillMissAllele = args[5])
