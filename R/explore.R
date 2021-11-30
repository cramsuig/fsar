## Removing stats_summary_table dependency -------------------------------------
source("https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/stats_summary_table.R")
mat <- mtcars; groups = setNames(rep(c("X1", "X2"), length.out = ncol(mat)), names(mat))
groups <- gtools::mixedsort(groups)
mat <- mat[, names(groups)]
stattab_sst <- stats_summary_table(
  mat = mat,
  groups = groups,
  moments = c('mn', 'sd')
)
stattab_red <- data.frame(cbind(
  t(apply(mat, 1, function(vec) tapply(vec, groups, matrixStats::mean2, na.rm = TRUE) )),
  t(apply(mat, 1, function(vec) tapply(vec, groups, stats::sd, na.rm = TRUE) ))
)); colnames(stattab_red) <- paste0(c("X1", "X2"), rep(c("_mean", "_sd"), each = 2))
all.equal(stattab_sst, stattab_red)

## fsar_gsea_test --------------------------------------------------------------
library(dplyr)
gse14308 <- GEOquery::getGEO("GSE14308")[[1]]
tmp <- setNames(gse14308@featureData$`Gene Symbol`, gse14308@featureData$ENTREZ_GENE_ID)
eranks <- setNames(fgsea::exampleRanks, casefold(tmp[names(fgsea::exampleRanks)], upper = TRUE))
msigdbr_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(grepl("HALLMARK", gs_name)) %>% split(x = .$gene_symbol, f = .$gs_name)
res_f <- fsar_gsea_test(eranks, msigdbr_sets, method = "fgsea")
res_f[, -ncol(res_f)]
res_l <- fsar_gsea_test(eranks, msigdbr_sets, method = "liger", verbose = TRUE)
res_l

# Resources --------------------------------------------------------------------
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
# https://bioconductor.org/packages/release/data/annotation/html/reactome.db.html
# https://www.biostars.org/p/375584/ # best ranking?
# Generate ranks:
# https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/

