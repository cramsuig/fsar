fsar_set_list = function(
  features,
  include = NULL,
  exclude = NULL,
  global = list(),
  fun = NULL, ...
) {
  features <- if(!is.list(features)) as.list(features) else features
  features <- lapply(features, FUN = function(x){
    if(file.exists(x[[1]][1]))
      x <- as.list(read.csv(x[[1]][1], stringsAsFactors = FALSE))
    if(is.list(x) || is.character(x)) x else as.list(x)
  })
  if(all(sapply(features, class) == "list"))
    features <- unlist(features, recursive = FALSE)
  features <- sapply(features, function(x){ y <- x[!is.na(x)]; y[which(y != "")] })
  features <- features[sapply(features, length) > 0]; #str(features)
  features <- c(global, features[!names(features) %in% names(global)])
  tmp = names(features)
  if(is.character(include)) tmp = tmp[grep(include, tmp, ignore.case = TRUE)]
  if(is.character(exclude)) tmp = tmp[!grepl(exclude, tmp, ignore.case = TRUE)]
  if(is.numeric(include)) tmp = tmp[include]
  if(is.numeric(exclude)) tmp = tmp[-exclude]
  features = features[tmp]
  if(!is.null(fun)) features = lapply(X = features, FUN = fun, ...)
  return(features)
}

fsar_scoring_set_lit <- function(
  features,
  grouping = FALSE,
  verbose = FALSE
){
  list_processed <- list()
  if(verbose) cat("Getting lists ready\n")
  for(i in 1:length(features)){
    mynameis <- names(features[i])
    tmp <- c(names(list_processed), sapply(list_processed, names))
    if(mynameis %in% tmp) next
    if(verbose) cat("@", mynameis, "\n")
    mylists <- features[i]
    if(isTRUE(grouping)){
      # First part of the name can be used to aggregate groups of list
      # e. g., tcell_treg, tcell_tfh, tcell_tfr
      tvar <- sub("^([[:alnum:]]{1,})_.*", "\\1", names(features[i]))
      tmp <- grep(tvar, names(features))
      if(length(tmp) > 1){
        if(verbose) cat(" - grouping in ", tvar, ": ", sep = "")
        mylists <- features[tmp]; mynameis <- tvar
        if(verbose) cat(names(mylists), "\n")
      }
    }#; names(mylists) <- paste0(names(mylists), ".score")
    list_processed[[mynameis]] <- c(list(name = mynameis), mylists)
  }; return(x = list_processed)
}
