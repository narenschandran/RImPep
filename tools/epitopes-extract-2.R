# This is used with a query FASTA with no TPM

library(getopt)
library(data.table)

spec <- matrix(c(
    'sample_dir'  , 'S', 1, 'character',
    'epitopes_dir', 'E', 1, 'character',


    'epitope_cutoff'     , 'e', 2, 'numeric',
    'epitope_rank_cutoff', 'r', 2, 'numeric',
    'presentation_cutoff', 'p', 2, 'numeric'
), byrow = T, ncol = 4)

opt <- getopt(spec)

stopifnot((!is.null(opt$epitopes_dir)) && (!is.null(opt$sample_dir)))
if (is.null(opt$epitope_cutoff)) opt$epitope_cutoff <- 0.7
if (is.null(opt$epitope_rank_cutoff)) opt$epitope_rank_cutoff <- 3
if (is.null(opt$presentation_cutoff)) opt$presentation_cutoff <- 0.9

sample_dir <- opt$sample_dir
al_f       <- file.path(sample_dir, 'alleles.txt')
pred_f     <- file.path(sample_dir, 'pred.tsv')
tpm_f      <- file.path(sample_dir, 'tpm.tsv')


alleles    <- readLines(al_f)

allele_to_str  <- function(x) gsub("[*:]", "-", x)

res_lst <- list()
for (allele in alleles) {
    allele_str <- allele_to_str(allele)
    if (!(allele_str %in% names(res_lst))) {
        epi_fname  <- paste0(allele_str, '.tsv.gz')
        epi_file   <- file.path(opt$epitopes_dir, epi_fname)

        epi_df     <- fread(epi_file, sep = '\t', header = T)
        best_epi_df <- subset(epi_df,
                              (bindprob >= opt$epitope_cutoff) &
                              (rank     <= opt$epitope_rank_cutoff))
        best_epi_df$allele <- allele
        res_lst[[allele]] <- best_epi_df
    }
}
res <- do.call("rbind.data.frame", res_lst)

write.table(res, file.path(sample_dir, 'epitopes.tsv'),
            sep = '\t', quote = F, row.names = F)
