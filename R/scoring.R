fsar_scoring <-  <- function(
  edata,
  metadata = NULL,
  prefix = NULL,
  signature_list = list(name = "Cell.Cycle", S = Seurat::cc.genes$s.genes, G2M = Seurat::cc.genes$g2m.genes),
  confounders = "RNA.*res",
  grouping = FALSE, # grouping lists into one NAME_description_source
  plotting = TRUE,
  verbose = FALSE
  ...
){
  if(verbose) cat("-- Signature scoring --\n")
  if(verbose) cat("Output:", prefix, "\n")
  edata <- if(casefold(edata) == "seurat"){
    Seurat::CreateSeuratObject(edata, meta.data = metadata)
  }
  if(!is.null(metadata)){
    edata@meta.data <- joindf(edata@meta.data, metadata)
  }

  if(any(!confounders %in% colnames(edata@meta.data))){
    confounders <- filters_columns(edata@meta.data, confounders, verbose = verbose)
  };
  confounders <- confounders[confounders %in% colnames(edata@meta.data)]
  if(verbose)
    cat("Confounders:", stringr::str_wrap(paste0(confounders, collapse = ", ")), "\n");
  Idents(edata) <- 'orig.ident'
  signature_class <- list()
  signature_class[["processed"]] <- feature_set_process(signature_list, grouping = grouping)
  signame <- paste0(prefix, "signatures.csv")
  if(file.exists(signame)){
    if(verbose) cat("Pre-computed scores\n")
    signdf <- readfile(signame, stringsAsFactor = FALSE, check.names = FALSE, row.names = 1)
    edata@meta.data <- joindf(edata@meta.data, signdf)
  }
  scores_cols <- unique(unlist(lapply(signature_class[["processed"]],
    function(x) c(x[[1]], names(x)[-1]) )))
  for(scoring in signature_class[["processed"]]){
    if(verbose) cat("\n@", scoring$name, "\n")
    check_both <- c(scoring$name, names(scoring[-1])) # check both names
    if(any(!check_both %in% colnames(edata@meta.data))){
      edata <- ClassifyScoring(
        edata = edata, name = scoring$name, verbose = verbose,
        features = scoring[head(2:length(scoring), 2)]
      ); str(scoring); str(scores_cols)
      tvar <- FetchData(edata, vars = scores_cols[scores_cols %in% colnames(edata@meta.data)])
      if(exists("signdf")) tvar <- joindf(tvar, signdf)
      write.csv(tvar, file = signame)
    }
  }
  signature_class[["scores"]] <- signdf
  if(isTRUE(plotting)) fsar_scoring_plot(signature_class, verbose = verbose, ...)
  if(verbose) cat("-- --------- ------- --\n")
  return(edata)
}
