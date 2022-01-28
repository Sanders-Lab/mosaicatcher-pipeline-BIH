# Whoeps, 07th Jan 2021
# Making a large overview over the results from the arbigent folder. 
# I'm giving myself 2h to make this nice today. 

# Input: callmatrix from clean_genotype.R
# Input: a csv from david from which to extract samplenames
# Output: a matrix with added entries:
#   Filter - Pass, NoReadsPass, MendelFail, FalsePositive, lowconf, AlwaysComplex
#   Mapability - percent?
#   nhom, nhet, nref, nnoreads, ncomplex
#   mendel 1/0

# Load libraries
library(stringr)
library(dplyr)
library(pheatmap)
library(matrixStats)
library(reshape2)
library(optparse)
source('postprocess_helpers.R')


#INPUT INSTRUCTIONS
option_list = list(
  make_option(c("-i", "--table"), type="character", default=NULL,
              help="res.csv", metavar="character"),
  make_option(c("-n", "--normal_names"), type="character", default=T,
              help="Should samplenames be converted from GM to NA?", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="Outfile: verdicted table", metavar="character")

)

# Parse input
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
####local testing
#opt$table='/home/hufsah/HILBERT_gpfs/HGSVC/ARBIGENT_2022/pipeline/regenotyper_allsamples_bulk/arbigent_results/res.csv'
#opt$normal_names<-1
#opt$outfile<-'/home/hufsah/HILBERT_gpfs/HGSVC/ARBIGENT_2022/pipeline/regenotyper_allsamples_bulk/arbigent_results/res_verdicted.vcf'
callmatrix_link = opt$table
normal_names = opt$normal_names
outfile = opt$outfile

#callmatrix_link = '~/s/g/korbel2/StrandSeq/Test_WH/pipeline_7may/pipeline/regenotyper_allsamples_bulk/arbigent_results/res.csv'
#normal_names = T

cm = read.table(callmatrix_link, header=1, sep='\t', stringsAsFactors = F )


### GO ###

# some inventory. Which samples do we have here? And therefore how many 'other' cols?
samples = colnames(cm)[grep("^[HMNG].*", colnames(cm))]
n_samples = length(samples)
n_other_cols = dim(cm)[2] - n_samples

# Rename samples if wanted
if (as.numeric(normal_names)){
  colnames(cm)[(n_other_cols+1):dim(cm)[2]] = 
    str_replace(colnames(cm)[(n_other_cols+1):dim(cm)[2]], 'GM', 'NA')
  samples = colnames(cm)[grep("^[HMNG].*", colnames(cm))]
}


# Factor char stuff
cm[] <- lapply(cm, as.character)

# stratify entries with 0 valid bins. 
cm[cm$valid_bins==0,samples] = 'noreads'

# Count hom, het, ref, noreads and complex
cm = count_homhetrefetc(cm, n_samples)

# Calc mapability
cm$valid_bins = as.numeric(cm$valid_bins)

# Mendel
#Hufsah 28.1.22
#Only run mendel Check if all trios are present
if ("NA19238" %in% colnames(cm)& "NA19239" %in% colnames(cm)&
    "NA19240" %in% colnames(cm)& "HG00512" %in% colnames(cm)&
    "HG00513" %in% colnames(cm)& "HG00514" %in% colnames(cm)&
    "HG00731" %in% colnames(cm) & "HG00732" %in% colnames(cm)&
    "HG00733" %in% colnames(cm)){
  cm = add_mendelfails(cm)
}else{
  cm$mendel1<-NA
  cm$mendel2<-NA
  cm$mendel3<-NA
  cm$mendelfails=NA
}

# Filter
if(length(samples)>1){
  cm = apply_filter_new(cm, samples)
}else{
  cm=apply_filter_onesample(cm,samples)
}
  

# Clean
cm[,c('mendel1','mendel2','mendel3')] = NULL
cm[] <- lapply(cm, as.character)

# Rename inv_dup genotypes
cm = make_invdups_human_readable(cm, n_samples)

# Sort columns
cols = c(colnames(cm)[1:n_other_cols],'verdict', 'nref','nhet','nhom','ninvdup','ncomplex',samples)
cm_return = cm[,cols]

# Save
write.table(cm_return, file=outfile, col.names=T, row.names=F, sep='\t', quote=F)



