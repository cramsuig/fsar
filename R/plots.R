fsar_scoring_plot <- function(
  results,
  confounders = "RNA.*res|^orig\\.|\\.tag$",
  reductions = list(umap = c('UMAP_1', 'UMAP_2')),
  scatter_col = NULL,
  violins_fun = if(exists("violin")){
    if(is.function(violin)) function(data, x, y) violin(data, x, y, colour_by = "pct")
  },
  return_plot = FALSE,
  output = "./",
  verbose = FALSE
) {
  if(is.null(scatter_col))
    scatter_col <- c("#fffffa", "#fffeee", "#ffe080", "#ffc100",
      "#ff0000", "#EE0000", "#a10000", "#670000")
  if(any(!confounders %in% colnames(results[["scores"]]))){
    confounders <- grep(
      pattern = confounders, x = colnames(results[["scores"]]), value = TRUE)
  };
  if(verbose)
    cat("Confounders:", stringr::str_wrap(paste0(confounders, collapse = ", ")), "\n");
  cols <- c(unname(unlist(reductions)), confounders)
  pp_list <- list()
  for(scoring in results[["processed"]]){
    ddfplot <- results[["scores"]][, c(scoring$name, names(scoring[-1]), cols)]
    ddfplot$Signature <- as.character(ddfplot[, scoring$name])
    ddfplot$Signature <- ifelse(ddfplot$Signature == "None", "No", "Yes")
    if(verbose) cat(" # heatmap\n")
    fname0 <- paste0(output, scoring$name)
    fname <- paste0(fname0, '_heatmap.pdf')
    if(!(file.exists(fname) && file.size(fname) > 3620)){
      tvar <- withr::with_seed(42,
        ddfplot[sample(1:nrow(ddfplot), min(c(1000, nrow(ddfplot)))), ])
      tvar <- tvar[order(tvar$Signature), ]
      p <- try(pheatmap::pheatmap(
        mat = data[unname(unlist(scoring[-1])), rownames(tvar)],
        scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
        annotation_col = tvar[, c(confounders, "Signature")],
        show_colnames = FALSE, silent = TRUE
      ))
      pdf(fname, width = 10, height = 12, onefile = FALSE); print(p); dev.off()
      pp_list[[fname]] <- p
    }
    for(i in names(reductions)){
      if(verbose) cat(" #", i, "\n")
      fname <- paste0(fname0, '_', i, '.pdf')
      if(file.exists(fname) && file.size(fname) > 3620) next
      p <- lapply(names(scoring[-1]), function(x){
        aesy <- aes_string(x = reductions[[i]][1], y = reductions[[i]][2], color = x)
        ggplot(data = ddfplot, mapping = aesy) +
          geom_point(size = 0.1) + scale_color_gradientn(colours = scatter_col) +
          labs(colour = NULL, title = x)
      })
      tvar <- c(10, 10)#plot_size(make_grid(length(p)))
      p <- cowplot::plot_grid(plotlist = p)
      pdf(fname, width = tvar[2], height = tvar[1]); print(p); graphics.off()
      pp_list[[fname]] <- p
    }
    if(is.null(violins_fun)) next
    for(confy in confounders){
      if(verbose) cat(" -", confy, "\n")
      fname <- paste0(fname0, "_violin_", confy, '.pdf')
      if(file.exists(fname) && file.size(fname) > 3620) next
      p <- violins_fun(ddfplot, confy, names(scoring[-1])) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      pdf(fname, width = 11, height = 7); print(p); graphics.off()
      pp_list[[fname]] <- p
    }
  }; if(isTRUE(return_plot)) return(pp_list) else invisible(x = NULL)
}

fsar_gsea_test_plot_gggsea <- function(
  results, rank, features, ...
) {
  features <- features[names(features) %in% results$pathway]
  results <- results[results$pathway %in% features_names, ]
  df <- gggsea::gseaCurve(rank, features[unique(results$pathway)], results)
  ggplot2::ggplot() + call_formal_args(gggsea::geom_gsea, list(df = df), ...)
}

fsar_gsea_test_plot <- function(
  results, rank, features,
  method = c("gggsea", "fgsea", "liger"),
  verbose = FALSE,
  output = "./",
  return_plot = TRUE
  ...
) {
  method <- match.arg(method)
  pp_list <- list()
  for (i in names(features)) {
    if(method == "gggsea"){
      p <- fsar_gsea_test_plot_gggsea(
        results, rank, features[i], ...
      )
    }
    # The rest I need to solve how to capture plot() into a grob
    fname <- paste0(output, i, ".pdf")
    pdf(fname); print(p); graphics.off()
    pp_list[[fname]] <- p
  }; if(isTRUE(return_plot)) return(pp_list) else invisible(x = NULL)
}
