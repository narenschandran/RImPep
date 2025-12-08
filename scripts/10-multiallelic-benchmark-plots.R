library(scriptloc)
source(file.path(script_dir_get(), '..', '.conf.R'))

library(ROCit)

almat_get <- function(datf) {
    al_pat  <- "^al[1-9]:"
    al_cols  <- grep(al_pat, colnames(datf), value = T)
    almat    <- if (length(al_cols) == 6) {
        almat   <- as.matrix(datf[, al_cols])
    } else {
        alleles <- sub(al_pat, "", al_cols)
        genes   <- substr(alleles, 1, 5)
        misgenes <- local({
            tbl <- table(genes)
            names(which(tbl == 1))
        })
        miscols <- al_cols[match(genes, misgenes)]
        as.matrix(datf[,c(al_cols, miscols)])
    }
}


base_d <- file.path(RES_DIR, '03-multiallele-benchmark')

ds <- local({
    tmp <- list.dirs(base_d, recursive = F)
    setNames(tmp, basename(tmp))
})

pred_fs <- setNames(file.path(ds, 'pred.tsv'),
                    names(ds))

stopifnot(all(file.exists(pred_fs)))
dat_lst <- lapply(pred_fs, function(pred_f) {
    x <- read.table(pred_f, sep = '\t', header = T, check.names = F,
                    stringsAsFactors = F)
    y <- local({
        tmp <- subset(x, label != 0.5)
        tmp$label[tmp$label > 0.5] <- 1
        tmp
    })
    list(datf = y, almat = almat_get(y))
})


res_df <- data.frame(
    sample = names(pred_fs)
)

proc_fn <- function(dat, fn) {
    y       <- dat$datf
    almat   <- dat$almat
    y$score <- fn(almat)
    with(y, rocit(score = score, class = label)$AUC)
}

res_lst <- list()
fn_lst <- list(
    'sum'   = rowSums,
    'max'   = function(x) apply(x, 1, max),
    'sqsum' = function(x) {
        rowSums(x * x)
    },
    'sqrt'  = function(x) {
        rowSums(sqrt(x))
    },
    'rank'  = function(x) {
        rnk   <- apply(-x, 2, rank, ties.method = 'max')
        rnkp <- (1 - apply(rnk, 2, function(y) y / max(y)))
        rowSums(rnkp)
    },
    'rank2'  = function(x) {
        rnk   <- apply(-x, 2, rank, ties.method = 'max')
        rnkp <- (1 - apply(rnk, 2, function(y) y / max(y)))
        rowSums(x * (rnkp > 0.9))
    },
    'rank3'  = function(x) {
        rnk   <- apply(-x, 2, rank, ties.method = 'max')
        rnkp <- (1 - apply(rnk, 2, function(y) y / max(y)))
        rowSums(x * (rnkp > 0.8))
    },
    'rank4'  = function(x) {
        rnk   <- apply(-x, 2, rank, ties.method = 'max')
        rnkp <- (1 - apply(rnk, 2, function(y) y / max(y)))
        rowSums(x * (rnk <= 2000))
    },
    'ranksub1'  = function(x) {
        rnk   <- apply(-x, 2, rank, ties.method = 'max')
        rnkp <- (1 - apply(rnk, 2, function(y) y / max(y)))
        rowSums(rnkp * (rnkp > 0.8))
    },
    'ranksub2'  = function(x) {
        rnk   <- apply(-x, 2, rank, ties.method = 'max')
        rnkp <- (1 - apply(rnk, 2, function(y) y / max(y)))
        rowSums(rnkp * (rnkp > 0.5))
    },
    'geom'  = function(x) {
        apply(x, 1, function(y) prod(y) ^ (1/length(y)))
    },
    'top3sum'  = function(x) {
        apply(x, 1, function(y) sum(sort(y, decreasing = T)[1:3]))
    },
    'top2sum'  = function(x) {
        apply(x, 1, function(y) sum(sort(y, decreasing = T)[1:2]))
    },
    'top4sum'  = function(x) {
        m <- t(apply(x, 1, sort, decreasing = T))
        rowSums(m[,1:4])
    }
    # 'top5sum'  = function(x) {
    #     apply(x, 1, function(y) sum(sort(y, decreasing = T)[1:5]))
    # }
)

system.time({
for  (fn_name in names(fn_lst)) {
    if (!(fn_name %in% names(res_lst))) {
        res_lst[[fn_name]] <- sapply(dat_lst, proc_fn, fn = fn_lst[[fn_name]])
    }
}
})
print(sapply(res_lst, quantile))


res_df$AUC <- res_lst[["sum"]]

bst_df <- subset(res_df, grepl("^bassani", sample))
map_f1 <- file.path(PREREQ_DIR, 'bassani-sternberg.tsv.gz')
map1   <- local({
    tmp <- read.table(map_f1, sep = '\t', header = T)
    setNames(unique(tmp[,c("sample_id", "sample")]),
             c("sample_id", "cell_line"))
})
bst_df$cell_line <- map1$cell_line[match(bst_df$sample, map1$sample_id)]

bst_v <- setNames(bst_df$AUC, bst_df$cell_line)

{
svg(file.path(PLOTS_DIR, 'bassani-sternberg-auc.svg'),
    height = 7, width = 4)
barplot(bst_v, ylim = c(0, 1), ylab = "AUC",
        col = hcl.colors(length(bst_v), "Harmonic"))
abline(h = 0.5, lty = 2)
dev.off()
}

hlt_df <- subset(res_df, grepl("^hlatlas", sample))

map_f2 <- file.path(PREREQ_DIR, 'hlatlas.map')
map2 <- read.table(map_f2, sep = '\t', header = T, stringsAsFactors = F)

hlt_inds <- match(hlt_df$sample, map2$sample_id)

capitalize <- function(x) sapply(strsplit(x, " "), function(y) {
    substr(y, 1, 1) <- toupper(substr(y, 1, 1))
    paste(y, collapse = " ")
})
hlt_df[,c("tissue", "donor")] <- map2[hlt_inds, c("tmap", "donor")]
hlt_df$tissue <- capitalize(hlt_df$tissue)

tsplit <- with(hlt_df, split(AUC, tissue))
{
svg(file.path(PLOTS_DIR, 'hlatals-tissue-auc.svg'))
par(mar = c(8, 3, 3, 3))
boxplot(tsplit, ylim = c(0.65, 1), 
        col = hcl.colors(length(tsplit), "Spectral"),
        las = 2)
dev.off()
}

hlt_mat <- local({
    tmp <- tapply(hlt_df$AUC, hlt_df[,c("tissue", "donor")], identity)
    tmp1 <- cbind(tmp, Median = apply(tmp, 1, median, na.rm = T))
    tmp2 <- rbind(tmp1,
                  Median = apply(tmp1, 2, median, na.rm = T))
    tmp2["Median", "Median"] <- median(tmp, na.rm = T)
    tmp2
})
hlt_ch  <- apply(hlt_mat, 2, function(v) {
    ch <- rep("", length(v))
    ch[!is.na(v)] <- sprintf("%0.2f", v[!is.na(v)])
    ch
})

library(pheatmap)

{
svg(file.path(PLOTS_DIR, 'hlatals-all-auc.svg'),
    height = 7.2, width = 7.2)
hlt_br  <- seq(0.75, 0.9, 0.005)
hlt_cl  <- rev(hcl.colors(length(hlt_br) - 1, "RdYlBu"))
hlt_ncl <- matrix("black", ncol = ncol(hlt_mat), nrow = nrow(hlt_mat))
hlt_ncl[hlt_mat > 0.87] <- "white"
pheatmap(
    hlt_mat,
    cluster_rows = F,
    cluster_cols = F,
    display_numbers = hlt_ch,
    na_col = "white",
    border_col = "black",
    col = hlt_cl,
    breaks = hlt_br,
    number_col = hlt_ncl,
    legend = F,
    cellwidth = 16,
    cellheight = 16,
    gaps_row = c(26, 26),
    gaps_col = c(21, 21)
)
dev.off()
}


{
svg(file.path(PLOTS_DIR, 'bassani-sternberg-alternative.svg'))
pheatmap(
    t(bst_v),
    cluster_rows = F,
    cluster_cols = F,
    display_numbers = T,
    na_col = "white",
    border_col = "black",
    col = hlt_cl,
    breaks = hlt_br,
    number_col = "black",
    legend = F,
    cellwidth = 16,
    cellheight = 16,
)
dev.off()
}
