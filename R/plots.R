fsar_scoring_plot <- function(
  results,
  reductions = list(pca = c('PC_1', 'PC_2'), tsne = c('tSNE_1', 'tSNE_2'), umap = c('UMAP_1', 'UMAP_2')),
  couls = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000"),
  violins_color = "mean",
  verbose = FALSE
) {
  for(scoring in results[["processed"]]){
    tmp <- c(scoring$name, names(scoring[-1]), unname(unlist(reductions)), confounders)
    ddfplot <- results[["scores"]][, tmp]
    ddfplot$Signature <- as.character(ddfplot[, scoring$name])
    ddfplot$Signature <- ifelse(ddfplot$Signature == "None", "No", "Yes")
    # for(i in names(scoring[-1])) ddfplot[, i] <- scales::rescale(ddfplot[, i], to = c(0, 1))
    tvar <- table(ddfplot[, 'Signature'])
    subtitl <- paste0(paste0(names(tvar), ": ", unname(tvar)), collapse = "; ")
    if(verbose) cat(" # heatmap\n")
    fname <- paste0(prefix, scoring$name, '_heatmap.pdf')
    if(!(file.exists(fname) && file.size(fname) > 3620)){
      pdf(fname, width = 10, height = 12, onefile = FALSE);
      p <- try(custom_heatmap(
        object = object,
        rnames = unname(unlist(scoring[-1])),
        orderby = scoring$name,
        sample_it = c("orig.ident", '3000'),
        categorical_col = c(confounders, "Signature"),
        feature_order = "pca",
        verbose = FALSE,
        show_colnames = FALSE
      ))
      graphics.off()
    }
    for(i in names(reductions)){
      if(verbose) cat(" #", i, "\n")
      fname <- paste0(prefix, scoring$name, '_', i, '.pdf')
      if(file.exists(fname) && file.size(fname) > 3620) next
      p <- lapply(names(scoring[-1]), function(x){
        aesy <- aes_string(x = reductions[[i]][1], y = reductions[[i]][2], color = x)
        ggplot(data = ddfplot, mapping = aesy) +
        geom_point(size = 0.1) + scale_color_gradientn(colours = couls) +
        labs(colour = NULL, title = make_title(x))
      })
      tvar <- plot_size(make_grid(length(p))) # t o determine when it's more
      pdf(fname, width = tvar[2], height = tvar[1]);
      print(cowplot::plot_grid(plotlist = p)); graphics.off()
    }
    for(confy in confounders){
      if(verbose) cat(" -", confy, "\n")
      fname <- paste0(
        prefix, scoring$name,
        "_violin_", violins_color, "_", confy, '.pdf'
      )
      if(file.exists(fname) && file.size(fname) > 3620) next
      p <- violins(
        dat = ddfplot,
        xax = confy,
        yax = names(scoring[-1]),
        colour_by = violins_color
      ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      pdf(fname, width = 12, height = 8);
      print(p)
      graphics.off()
    }
  }
}
