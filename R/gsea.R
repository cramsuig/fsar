#' Ranking Metrics
#'
#' @description `fsar_gsea_metric` creates values useful for ranking in a GSEA.
#'
#' @details This is a function that takes a matrix of measurements and uses it to
#' calculate metric values from mean, sd, or a test that will be used by `gsea_tests`.
#'
#' @param mat Count matrix.
#' @param groups Column name containing the groups.
#' @param metric Type of metric.
#' @param rnames Feature names to use.
#' @param verbose Show progress.
#' @keywords mean sd
#'
#' @return A vector of values.
#'
#' @author Ciro Ramírez-Suástegui, \email{ksuasteguic@gmail.com}
#' @references \url{https://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html}
#' @seealso \code{\link{gsea_tests}}
#'
#' @importFrom stringr str_wrap
#' @importFrom gtools mixedsort
#' @importFrom matrixStats mean2 rowMaxs
#' @importFrom stats sd
#' @importFrom base rowSums
#'
#' @examples
#' gsea_ready = fsar_gsea_metric(mat = edata, groups = c(rep("G1", 200), rep("G2", 200)))
#'
#' @export
#'

fsar_gsea_metric <- function(
  mat,
  groups,
  metric = c('Signal2Noise', 'tTest', 'Ratio_of_Classes', 'Diff_of_Classes', 'log2_Ratio_of_Classes'),
  verbose = FALSE
) {
  metric <- match.arg(metric)
  tvar <- unique(groups) # first to appear will be "up"
  if(verbose) cat("Groups:", paste0(tvar, collapse = ", "), "\n")
  if(verbose) cat("X2 =", tvar[2], "\n")
  groups <- sub(paste0("^", tvar[1], "$"), "X1", groups)
  groups <- sub(paste0("^", tvar[2], "$"), "X2", groups)
  groups <- gtools::mixedsort(groups)
  mat <- mat[, names(groups)]
  stattab <- data.frame(cbind(
    t(apply(mat, 1, function(vec) tapply(vec, groups, matrixStats::mean2, na.rm = TRUE) )),
    t(apply(mat, 1, function(vec) tapply(vec, groups, stats::sd, na.rm = TRUE) ))
  )); colnames(stattab) <- paste0(c("X1", "X2"), rep(c("_mean", "_sd"), each = 2))
  stattab2 <- data.frame(sapply(colnames(stattab), function(x){
    y <- stattab[, x]
    # σ has a minimum value of .2 * absolute(μ), where μ=0 is adjusted to μ=1
    if(grepl("_sd", x)){
      mean_adj = abs(stattab[, sub("_sd", "_mean", x)])
      mean_adj = ifelse(mean_adj == 0, 1, mean_adj)
      matrixStats::rowMaxs( as.matrix(data.frame(minsd = .2 * mean_adj, sd = y)) )
    }
    return(y)
  }), row.names = rownames(stattab)); stattab <- stattab2; rm(stattab2)
  means <- rev(grep("_mean", colnames(stattab)))
  if(verbose) cat("Means:", means, "\n")
  dmean <- apply( stattab[, means], 1, diff )
  if(verbose) cat(metric, "\n")
  ymetric <- switch(metric,
    Signal2Noise = {
      if(verbose) cat("X1mean - X2mean\n---------------\n  X1sd + X2sd\n")
      dmean / apply( stattab[, grep("_sd", colnames(stattab))], 1, sum )
    }, tTest = {
      if(verbose) cat("          X1mean - X2mean\n")
      if(verbose) cat("  ------------------------------\n")
      if(verbose) cat("sqrt ((X1sd^2 / X1n) +( X2sd^2 / X2n))\n")
      ttab <- sweep(stattab[, grep("_sd", colnames(stattab))] ** 2, 2, table(groups), "/")
      dmean / sqrt( base::rowSums(ttab) )
    }, Ratio_of_Classes = {
      if(verbose) cat("X1mean / X2mean\n")
      stattab[, "X1_mean"] / stattab[, "X2_mean"]
    }, Diff_of_Classes = {
      if(verbose) cat("X1mean - X2mean\n"); dmean
    }, log2_Ratio_of_Classes = {
      if(verbose) cat("log2(X1mean / X2mean)\n")
      log2(stattab[, "X1_mean"] / stattab[, "X2_mean"])
    }
  )
  ymetric
}

fsar_gsea_test <- function(
  rank,
  features,
  method = c("fgsea", "liger"),
  verbose = FALSE,
  ...
) {
  method <- match.arg(method)
  if(verbose) cat("Method:", method, "\n")
  if(verbose) cat("Sets:", length(features), "\n")
  if(method == "liger"){
    results <- liger::bulk.gsea(values = rank, set.list = features, ...)
    colnames(results) <- c("pval", "padj", "ES", "edge")
    results <- cbind(pathway = rownames(results), results)
    results$size = sapply(features, function(x) sum(x %in% names(rank)) )
  }
  if(method == "fgsea"){
    results <- fgsea::fgsea(pathways = features, stats = rank, ...)
    # func <- fgsea::fgsea
    # options <- list(pathways = features, stats = rank)
    # options[["nproc"]] = min(c(ceiling(parallel::detectCores() / 5), 4))
  }; results$method = method
  tibble::as_tibble(results)
}

#' @title Call a function with its formal arguments.
#' @description Call a function adapting the formal arguments
#' @param fun Function.
#' @param opt_pre List of arguments to modify.
#' @param ... Extra arguments.
#' @return 'fun' output.
#' @details We first append ... with 'opt_pre' and then match it with the formal
#' arguments present in 'fun'.
#' @examples
#' \dontrun{
#' if(interactive())
#'   result <- call_formal_args(mean, list(x = rnorm(5), na.rm = FALSE))
#' @seealso
#'  \code{\link[methods]{methodUtilities}}
#'  \code{\link[withr]{with_seed}}
#' @rdname call_formal_args
#' @export
#' @importFrom methods formalArgs
#' @importFrom withr with_seed

call_formal_args <- function(fun, opt_pre, ...){
  opt_glo <- list(...)
  opt <- c(opt_pre, opt_glo[setdiff(names(opt_glo), names(opt_pre))])
  opt <- opt[intersect(names(opt), methods::formalArgs(fun))]
  withr::with_seed(42, { do.call(what = fun, args = opt) })
}
