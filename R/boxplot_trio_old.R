
#' @importFrom grDevices col2rgb
#' @importFrom stats setNames

boxplot_trio_v1 <- function(polar,
                         value,
                         box_colours = c('green3', 'blue', 'red'),
                         test = "polar_pvalue",
                         levels_order = NULL,
                         my_comparisons = NULL,
                         text_size = 10,
                         stat_colour = "black",
                         stat_size = 3,
                         step_increase = 0.05,
                         plot_method="ggplot",
                         ...){

  sampledata <- polar@sampledata
  expression <- polar@expression
  pvalues <- polar@pvalues

  if(! test %in% c("polar_pvalue", "polar_padj", "polar_multi_pvalue",
                   "polar_multi_padj", "t.test", "wilcox.test", "anova",
                   "kruskal.test")) {
    stop(paste("expression must be a data frame or c('polar_pvalues',",
               "'polar_padj', 'polar_multi_pvalue', 'polar_multi_padj',",
               "'t.test', 'wilcox.test', 'anova', 'kruskal.test')"))
  }
  if(is.null(levels_order)) {
    levels_order <- levels(sampledata[, polar@contrast])
  }
  if(! inherits(levels_order, "character")) {
    stop("levels_order must be a character vector")
  }
  if(! all(levels_order %in% levels(sampledata[, polar@contrast]))){
    stop(paste('levels_order must be a character vector defining the order',
               'of levels in sampledata[, contrast]'))
  }
  if(length(box_colours) != length(levels_order)){
    stop(paste0('The length of box_colours must match teh length of',
                'levels_order'))
  }
  # Convert to hex colours for plotly
  box_colours <- unlist(lapply(box_colours, function(x) {
    if(! grepl("#", x) &
       class(try(col2rgb(x), silent = TRUE))[1] == "try-error") {
      stop(paste(x, 'is not a valid colour'))
    } else if (! grepl("#", x) ) {
      y <- col2rgb(x)[, 1]
      x <- rgb(y[1], y[2], y[3], maxColorValue=255)
    }
    return(x)
  }))

  if(! is.data.frame(sampledata)) {
    stop("sampledata must be a data frame")
  }
  if(! inherits(expression, c("data.frame", "matrix"))) {
    stop("expression must be a data frame or matrix")
  }
  if(! inherits(value, c("character", "numeric"))) {
    stop("value must be a character")
  }
  if(length(value) > 1) stop("value must be of length 1")
  if(! value %in% rownames(expression)) {
    stop("value/gene is not in rownames(expression)")
  }
  if(! identical(colnames(expression), as.character(sampledata$ID))) {
    stop("expression and sampledata misalligned")
  }
  if(grepl("multi", test) & is.null(polar@multi_group_test)){
    stop(paste("A multi-group test parameter is required in pvalues to use",
               test))
  }
  if(! plot_method %in% c('plotly', 'ggplot')){
    stop("plot_method must be either plotly or ggplot")
  }

  colour_map <- setNames(box_colours, levels_order)
  sampledata$comp <- sampledata[, polar@contrast]
  expression <- expression[, match(as.character(sampledata$ID),
                                   colnames(expression))]

  if(is.character(value)) {
    index <- which(rownames(expression) ==  value)
  }

  if(is.null(my_comparisons)) {
    comps <- levels_order
    my_comparisons <- lapply(seq_len(ncol(combn(comps, 2))), function(i) {
      as.character(combn(comps, 2)[, i])
    })
  }

  df <- data.frame("ID" = sampledata$ID,
                   "group" = sampledata$comp,
                   "row" = as.numeric(as.character(expression[value, ])))
  df <- df[! is.na(df$row), ]
  df <- df[df$group %in% levels_order, ]

  # relevel based on defined order
  df$group <- factor(df$group, levels_order)
  df$col <- factor(df$group, labels=colour_map[match(levels(df$group),
                                                     names(colour_map))])

  map_pos <- setNames(seq_along(levels(df$group)), levels(df$group))

  # Calculate the pvalues depending on test
  if(test %in% c("t.test", "wilcox.test", "anova", "kruskal.test")){
    pvals <- compare_means(formula = row ~ group, data = df,
                           comparisons = my_comparisons,
                           method = test,
                           step.increase = step_increase,
                           size=stat_size)

    pvals$x.position <- map_pos[pvals$group1] +
      (map_pos[pvals$group2] - map_pos[pvals$group1])/2
    pvals$y.position <- max(df$row, na.rm=TRUE)*
      (1.01 + step_increase*c(seq_len(nrow(pvals))-1))
    pvals$new_p_label <- pvals$p.format

    # groups comparisons
  } else if (! grepl("multi", test)){
    if(! any(grepl(gsub("polar_", "", test), colnames(pvalues)))){
      stop(paste(test, "tests must have", gsub(".*_", "", test), 
                 "columns in polar@pvalues"))
    }
    pvals <- pvalues[value, ]
    pvals <- pvals[, grepl(gsub("polar_", "", test), colnames(pvals))]
    if(! is.null(polar@multi_group_test)){
      pvals <- pvals[, ! grepl(polar@multi_group_test, colnames(pvals))]
    }
    colnames(pvals) <- gsub(paste0("_", gsub("polar_", "", test)), "",
                                 colnames(pvals))
    rownames(pvals) <- "p"
    pvals <- data.frame(t(pvals))
    pvals$group1 <- gsub("_.*", "", rownames(pvals))
    pvals$group2 <- gsub(".*_", "", rownames(pvals))
    pvals$p.format <- format(pvals$p, digits=2)
    pvals$method <- test
    pvals$y.position <- max(df$row, na.rm=TRUE)
    pvals$comp <- paste0(pvals$group1, "_", pvals$group2)
    pvals$x.position <- map_pos[pvals$group1] +
      (map_pos[pvals$group2] - map_pos[pvals$group1])/2
    pvals$y.position <- pvals$y.position[1]*
      (1.01 + step_increase*(seq_len(nrow(pvals))-1))
    comp_use <- unlist(lapply(my_comparisons, function(x) { 
      c(paste0(x[1], "_", x[2]), paste0(x[2], "_", x[1]))
      }))
    pvals <- pvals[pvals$comp %in% comp_use, ]

    # muti group comparisons
  } else{
    if(! any(grepl(gsub("polar_multi_", "", test), colnames(pvalues)))){
      stop(paste(test, "tests must have", gsub(".*_", "", test), 
                 "columns in polar@pvalues"))
    }
    pvals <- pvalues[value, ]
    pvals <- pvals[, grepl(gsub("polar_multi_", "", test), colnames(pvals))]
    pvals <- pvals[, grepl(polar@multi_group_test, colnames(pvals))]
  }

  if(plot_method == 'ggplot'){
    p <- ggboxplot(data = df,
                   x = "group",
                   y = "row",
                   xlab = "",
                   ylab = paste(polar@pvalues$label[index], "Expression"),
                   fill = "group",
                   color = "group",
                   palette = box_colours,
                   outlier.shape = NA,
                   alpha = 0.3) +
      geom_jitter(data=df, height = 0, width = 0.30,
                  aes_string(color="group")) +
      theme(legend.position = "none",
            text = element_text(size = text_size),
            plot.background = element_rect(fill="transparent", color=NA),
            panel.background = element_rect(fill="transparent", colour=NA),
            legend.background = element_rect(fill="transparent", colour=NA))


    if(! grepl("multi", test)){

      p <- p + stat_pvalue_manual(
        data = pvals, label = "p.format",
        xmin = "group1", xmax = "group2",
        step.increase = step_increase,
        y.position = "y.position", color = stat_colour,
        size=stat_size, ...)
    } else{
      # muti group comparison
      p <- p + annotate("text", x = 0.5 + length(unique(df$group))/2,
                        y = Inf, vjust = 2, hjust = 0.5, color = stat_colour,
                        label = paste("p =", format(pvals, digits = 2)))
    }
  } else{
    p <- df %>%
      plot_ly() %>%
      add_trace(x = ~as.numeric(group),  y = ~row,
                type = "box",
                colors = levels(df$col), color = ~col,
                opacity=0.5, marker = list(opacity = 0),
                hoverinfo="none", showlegend = FALSE) %>%
      add_markers(x = ~jitter(as.numeric(group)), y = ~row,
                  marker = list(size = 6, color=~col),
                  hoverinfo = "text",
                  text = ~paste0(ID,
                                 "<br>Group: ", group,
                                 "<br>Expression: ", row),
                  showlegend = FALSE) %>%
      layout(legend = list(orientation = "h",
                           x =0.5, xanchor = "center",
                           y = 1, yanchor = "bottom"
      ),
      xaxis = list(title = polar@contrast, tickvals = 1:3,
                   ticktext = levels(df$group)),
      yaxis = list(title = paste(value, "Expression")))


    lines <- list()
    if(! grepl("multi", test)){
      for (i in seq_len(nrow(pvals))) {
        line <- list(line=list(color = stat_colour))
        line[["x0"]] <- map_pos[pvals$group1][i]
        line[["x1"]] <- map_pos[pvals$group2][i]
        line[c("y0", "y1")] <- pvals$y.position[i]
        lines <- c(lines, list(line))
      }

      a <- list(
        x = as.numeric(pvals$x.position),
        y = pvals$y.position,
        text = format(pvals$p.format, digits=3),
        xref = "x",
        yref = "y",
        yanchor = "bottom",
        font = list(color = stat_colour),
        showarrow = FALSE
      )

      p <- p %>% layout(annotations = a, shapes=lines)

    } else{
      a <- list(
        x = 2,
        y = step_increase + max(df$row, na.rm=TRUE),
        text = format(pvals, digits=3),
        xref = "x",
        yref = "y",
        font = list(color = stat_colour),
        showarrow = FALSE
      )

      p <- p %>% layout(annotations = a)
    }


  }

  return(p)
}

