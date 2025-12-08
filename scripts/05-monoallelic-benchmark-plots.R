library(scriptloc)
source(file.path(script_dir_get(), '..', '.conf.R'))



base_d <- file.path(RES_DIR, '01-ablation-study', 'p3000')


ds <- local({
    tmp0 <- file.path(list.dirs(base_d, recursive = F), 'seed-00')
    tmp1 <- Filter(dir.exists, tmp0)
    setNames(tmp1, basename(dirname(tmp1)))
})


dat_lst <- lapply(ds, function(d) {
    mname <- basename(dirname(d))
    mseed <- basename(d)
    pred_f <- file.path(d, 'all-pred.tsv')
    conf_f <- file.path(d, 'conf.dct')

    conf <- local({
        x <- read.table(conf_f, sep = '\t', header = F)
        y <- setNames(as.data.frame(t(x[,2])), x[,1])
        y$model <- mname
        y$seed  <- mseed
        y
    })

    pred <- subset(read.table(pred_f, sep = '\t', header = T),
                   split_name == 'test')
    list(pred = pred, conf = conf)
})

conf_df <- do.call("rbind.data.frame", lapply(dat_lst, `[[`, "conf"))
rownames(conf_df) <- NULL

pats <- c(
    "[-]tse[-]" = "Full Model",
    "[-]tsx[-]" = "TPM + Seq",
    "[-]xse[-]" = "Seq + Epi",
    "[-]txe[-]" = "TPM + Epi",
    "[-]txx[-]" = "TPM",
    "[-]xsx[-]" = "Seq",
    "[-]xxe[-]" = "Epi"
)

conf_df$Name <- NA
for (pat in names(pats)) {
    conf_df$Name[grepl(pat, conf_df$model)] <- pats[[pat]]
}


library(ROCit)



# pats <- c(
#     "Full Model" = "",
#     "TPM + Seq"  = "",
#     "Seq + Epi"  = "",
#     "TPM + Epi"  = "",
#     "TPM"        = "",
#     "Seq"        = "",
#     "Epi"        = ""
# )

cls <- c(
    "Full Model" = "black",
    "TPM + Seq"  = "orange",
    "Seq + Epi"  = "blue",
    "TPM + Epi"  = "red",
    "TPM"        = "red",
    "Seq"        = "turquoise",
    "Epi"        = "purple"
)

ltys <- c(
    "Full Model" = "solid",
    "TPM + Seq"  = "dashed",
    "Seq + Epi"  = "solid",
    "TPM + Epi"  = "solid",
    "TPM"        = "dotted",
    "Seq"        = "dotdash",
    "Epi"        = "dotted"
)

if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = T)
{
svg(file.path(PLOTS_DIR, 'ablation-study.svg'))
par(cex = 1.5)
auc_vals <- c()
for (i in seq_along(pats)) {
    mod_nm <- pats[[i]]
    mod_id <- conf_df$model[conf_df$Name == mod_nm]
    stopifnot(length(mod_id) == 1)
    pred <- dat_lst[[mod_id]]$pred
    roc_empirical <- rocit(score = pred$pred, class = pred$label)
    if (i == 1) {
        plot(roc_empirical, col = c(1, cls[mod_nm]),
             legend = F, YIndex = F, lwd = 3, lty = ltys[[mod_nm]],
             xlim = c(0, 1), ylim = c(0, 1))
    } else {
        lines(roc_empirical$FPR, roc_empirical$TPR,
              col = cls[mod_nm], lwd = 3, lty = ltys[[mod_nm]])
    }
    auc_vals <- c(auc_vals, roc_empirical$AUC)
}
auc_vals <- setNames(auc_vals, pats)

legend("bottomright",
    sprintf("%s (AUC: %0.3f)", names(auc_vals), auc_vals),
    lty = ltys[names(auc_vals)],
    col = cls[names(auc_vals)], bty = 'n', lwd = 2, cex = 0.8)
dev.off()

}

fpred <- dat_lst[['model-tse-adcca349']]$pred

fpred_auc <- local({
    tmp <- sapply(split(fpred, fpred$allele), function(pred) {
        rocit(score = pred$pred, class = pred$label)$AUC
    })
    data.frame(
        allele = names(tmp),
        auc    = tmp
    )
})

{
svg(file.path(PLOTS_DIR, 'ablation-study-allelewise-auc.svg'), height = 5, width = 7)
par(cex = 1.5)
hist(
    fpred_auc$auc,
    xlab = "AUC",
)
dev.off()
}

s <- local({
    tmp <- dat_lst[["model-txx-0c57412b"]]$pred[,c("tpm", "pred")]
    tmp$pred <- round(tmp$pred, 3)
    unique(tmp)
})
s <- s[order(s$pred),]

plot(log1p(s$tpm), s$pred)

v <- local({
    tmp <- dat_lst[["model-xsx-21f3f3c6"]]$pred[,c("gene", "pred")]
    tmp$pred <- round(tmp$pred, 3)
    unique(tmp)
})
v <- v[order(v$pred),]


