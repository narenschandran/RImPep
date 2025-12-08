library(scriptloc)
script_dir <- script_dir_get()
projroot   <- file.path(script_dir, '..')
prereq_dir <- file.path(projroot, 'prereq')
data_dir   <- file.path(projroot, 'data')
res_dir    <- file.path(projroot, 'results')

fs <- file.path(prereq_dir, 
                c('bassani-sternberg.tsv.gz',
                  'hlatlas.tsv.gz'))

dat <- do.call("rbind",
               lapply(fs, read.table, sep = '\t', header = T))

dat_sp <- split(dat, dat$sample_id)

base_odir <- file.path(res_dir, '03-multiallele-benchmark')

for (nm in names(dat_sp)) {
    print(nm)
    sp <- dat_sp[[nm]]
    alleles_str <- unique(sp$alleles)
    stopifnot(length(alleles_str) == 1)

    alleles <- strsplit(alleles_str, ";")[[1]]
    tpm     <- sp[,c("gene", "tpm", "label")]

    odir <- file.path(base_odir, nm)
    if (!dir.exists(odir)) dir.create(odir, recursive = T)
    write.table(tpm, file.path(odir, 'tpm.tsv'),
                sep = '\t', row.names = F, quote = F)
    writeLines(alleles, file.path(odir, 'alleles.txt'))
}
