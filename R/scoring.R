fsar_scoring_classify <- function(
  object,
  features,
  set.ident = FALSE,
  name = 'Set',
  verbose = FALSE,
  ...
) {
  cc.columns <- grep(
    pattern = name, x = colnames(x = object@meta.data), value = TRUE)
  if(is.null(names(features))) names(features) <- paste0("S", 1:length(features))
  if(name %in% cc.columns){ warning(name, " pre-computed"); return(object) }
  if(verbose) str(features)

  # Renaming columns colliding with previous signatures; first the 'name'
  cc.columns <- make.names(c(cc.columns, name), unique = TRUE);
  name <- cc.columns[length(cc.columns)]
  cc.columns <- make.names( # now the 'classes'
    c(colnames(x = object@meta.data), names(features)), unique = TRUE)
  names(features) <- cc.columns[tail(1:length(cc.columns), length(names(features)))]

  classes <- names(features)
  if(verbose) cat("Calculating scores\n")
  ctrl_n <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  if(ctrl_n > 1000) ctrl_n <- 100
  if(verbose) cat("Background size:", ctrl_n, "\n")
  object <- Seurat::AddModuleScore(
    object = object,
    features = features,
    name = name,
    ctrl = ctrl_n,
    ...
  )
  cc.columns <- grep(
    pattern = name, x = colnames(x = object@meta.data), value = TRUE)
  cc.scores <- object@meta.data[, cc.columns]
  object@meta.data <- object@meta.data[, !colnames(object@meta.data) %in% cc.columns]
  if(verbose) cat("Classification based on score.\nNone: all < 0; Undecided: max score > 1.\n")
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores) {
      if (all(scores < 0)) {
        return("None")
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Undecided')
        } else {
          return(classes[which(x = scores == max(scores))])
        }
      }
    }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c('rownames', classes, name)
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c(classes, name)]
  if(verbose) cat("Adding to meta data:", paste0(c(classes, name), collapse = ", "), "\n")
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    if(verbose) cat("Setting classes as identities\n")
    object[['old.ident']] <- Idents(object = object)
    Idents(object = object) <- 'Class'
  }
  return(object)
}

fsar_scoring <-  function(
  edata,
  metadata = NULL,
  output = "./",
  signature_list = list(name = "Cell.Cycle", S = Seurat::cc.genes$s.genes, G2M = Seurat::cc.genes$g2m.genes),
  grouping = FALSE, # grouping lists into one NAME_description_source
  plotting = FALSE,
  verbose = FALSE
  ...
){
  if(verbose) cat("---------- Signature scoring ----------\n")
  if(verbose) cat("Output:", output, "\n")
  edata <- if(casefold(edata) == "seurat"){
    Seurat::CreateSeuratObject(edata, meta.data = metadata)
  }
  if(!is.null(metadata)){
    edata@meta.data <- tibble::column_to_rownames(dplyr::left_join(
      tibble::rownames_to_column(edata@meta.data),
      tibble::rownames_to_column(metadata)))
  }
  Idents(edata) <- 'orig.ident'
  signature_class <- list()
  signature_class[["processed"]] <- feature_set_process(signature_list, grouping = grouping)
  signame <- paste0(output, "signatures.csv")
  if(file.exists(signame)){
    if(verbose) cat("Pre-computed scores\n")
    signdf <- read.csv(signame, row.names = 1)
    edata@meta.data <- tibble::column_to_rownames(dplyr::left_join(
      tibble::rownames_to_column(edata@meta.data),
      tibble::rownames_to_column(signdf)))
  }
  scores_cols <- unique(unlist(lapply(signature_class[["processed"]],
    function(x) c(x[[1]], names(x)[-1]) )))
  for(scoring in signature_class[["processed"]]){
    if(verbose) cat("\n@", scoring$name, "\n")
    check_names <- c(scoring$name, names(scoring[-1])) # check all names
    if(any(!check_names %in% colnames(edata@meta.data))){
      # classify is actually only useful when there are several lists in one name
      edata <- fsar_scoring_classify(
        edata = edata, name = scoring$name, verbose = verbose,
        features = scoring[head(2:length(scoring), 2)]
      ); tvar <- intersect(scores_cols, colnames(edata@meta.data))
      if(exists("signdf"))
        tvar <- tibble::column_to_rownames(dplyr::left_join(
          tibble::rownames_to_column(edata@meta.data[, tvar]),
          tibble::rownames_to_column(signdf)))
      write.csv(tvar, file = signame)
    }
  }
  signature_class[["scores"]] <- signdf
  if(isTRUE(plotting))
    fsar_scoring_plot(signature_class, output = output, verbose = verbose, ...)
  if(verbose) cat("---------- --------- ------- ----------\n")
  invisible(x = NULL)
}
