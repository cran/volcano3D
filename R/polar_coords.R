setClassUnion("character_or_NULL", c("character", "NULL"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))

#' An S4 class to define the polar coordinates and pvalues for polar
#' differential expression plots
#' 
#' @slot sampledata Sample data with ID and contrast column.
#' @slot contrast The column name in `sampledata`` which contains the
#'   three-group contrast factor used for comparisons.
#' @slot pvalues A data frame containing the p-values, and adjusted p-values,
#'   for all three comparisons between groups in the 
#'   contrast factor, as well as optional fold changes and multi-group tests.
#' @slot multi_group_test Column name prefix for statistical tests between
#'   all three groups
#' @slot expression A data frame or matrix containing the expression data. This 
#' is used to calculate z-score and fold change, therefore it should be a 
#' normalised expression object such as log transformed or variance stabilised. 
#' @slot polar A data frame containing:
#'   \itemize{
#'       \item The axis score or mean expression for each of the three groups 
#'       in comparison
#'       \item The z-score polar coordinates: 'y_zscore', 'x_zscore' and 
#'       'r_zscore'
#'       \item The fold-change polar coordinates: 'y_fc', 'x_fc' and 'r_fc'
#'       \item 'angle': The angle in radians for polar coordinates
#'       \item 'angle_degrees': The angle in degrees
#'       \item 'max_exp': The group with the highest expression
#'       \item 'sig': The significance category
#'   }
#' @slot non_sig_name The category name for variables which are not significant

setClass("polar", slots = list(sampledata = "data.frame",
                               contrast = "character",
                               pvalues = "data.frame",
                               multi_group_test = "character_or_NULL",
                               expression = "df_or_matrix",
                               polar = "df_or_matrix",
                               non_sig_name = "character"))


#' Coordinates for Three Way Polar Plot
#'
#' This function creates a polar object of S4 class for downstream plots 
#' containing the p-values from a three-way group comparison, expression data 
#' sample data and polar coordinates.
#' @param sampledata A data frame containing the sample information.
#' This must contain: an ID column containing the sample IDs which can be 
#' matched to the `expression` data and a 
#' contrast column containing the three-level factor used for contrasts.
#' @param contrast The column name in `sampledata` which contains the 
#' three-level factor used for contrast.
#' @param pvalues A data frame containing: \itemize{
#' \item three `p_col_suffix` columns: one for 
#' the pvalue for each comparison between groups;
#' \item three optional `fc_col_suffix` columns for the fold change between 
#' each comparison (if NULL, no Fold Change columns are included); 
#' \item three optional `padj_col_suffix` columns (if NULL 
#' adjusted p values are calculated using `padjust_method`); 
#' \item and optional 'p', 
#' 'padj and 'fc' columns for a three-way test, such as ANOVA or likelihood 
#' ratio test, defined by `multi_group_prefix`. 
#' }
#' @param expression An optional data frame containing expression data for
#' downstream analysis and visualisation. The rows must contain probes which
#' match the rows in pvalues and the columns must contain samples which match
#' \code{sampledata$ID}.
#' @param groups The groups to be compared (if NULL this defaults
#' to \code{levels(sampledata[, 'contrasts'])}).
#' @param p_col_suffix The suffix word to define columns containing p values
#' (default = 'pvalues'). These must not contain underscores; 
#' @param padj_col_suffix The suffix word to define columns containing adjusted
#' p values (default = 'padj').  These must not contain underscores. 
#' If NULL these will be calculated using \code{padjust_method} 
#' @param fc_col_suffix The optional suffix word to define columns containing 
#' log fold change values (default = 'logFC'). These must not contain 
#' underscores.
#' @param padjust_method The method used to calculate adjusted p values if 
#' padj_col_suffix is NULL (default = 'BH'). See \code{\link[stats]{p.adjust}}. 
#' If NULL no adjusted pvalue is calculated. 
#' @param multi_group_prefix Optional column prefix for statistics (p, padj, 
#' and fold change) across all three groups (typically ANOVA or likelihood 
#' ratio tests). default = NULL. These must not contain underscores
#' @param non_sig_name Category name to assign to non-significant points
#' @param significance_cutoff Value defining the significance cut-off (points 
#' with pvalues below this point will be classed as \code{non_sig_name})
#' @param fc_cutoff The cut-off for fold change, below which markers will be
#' classed as \code{non_sig_name}` (default = 0.3).
#' @param label_column Optional column name in pvalues for markers to be 
#' labelled with at plotting stage. If NULL the rownames of pvalues are used. 
#' @param cutoff_criteria Whether to use pvalue or padj for the colour coding 
#' significance cutoff. 
#' @return Returns an S4 polar object containing:
#' \itemize{
#'   \item{'polar'} A data.frame containing:
#'   \itemize{
#'       \item{The mean expression for each of the three groups in comparison}
#'       \item{The z-score polar coordinates: 'y_zscore', 'x_zscore' and 
#'       'r_zscore'}
#'       \item{The fold-change polar coordinates: 'y_fc', 'x_fc' and 'r_fc'}
#'       \item{'angle': The angle in radians for polar coordinates}
#'       \item{'angle_degrees': The angle in degrees}
#'       \item{'max_exp': The maximally expressed group}
#'       \item{'sig': The significance group}
#'   }
#'   \item{'pvalues'} A data frame containing the p-values, adjusted p-values,
#'   and optional log(fold changes) for all three comparisons between groups in 
#'   the contrast factor, as well as optional multi-group tests.
#'   \item{'sampledata'} Sample data with column ID and contrast
#'   \item{'contrast'} The column name in `sampledata`` which contains the
#'   three-group contrast factor used for comparisons.
#'   \item{'multi_group_test'} Column name prefix for statistical tests between
#'   all three groups
#'   \item{'expression'} An optional data frame or matrix containing the
#'   expression data
#'   \item{'non_sig_name'} The category name for variables which are classed as 
#'   not significant
#' }
#' @keywords dplot spatial
#' @importFrom grDevices col2rgb hsv
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @export
#' @examples
#' data(example_data)
#' syn_polar <- polar_coords(sampledata=syn_example_meta,
#'                     contrast="Pathotype",
#'                     groups = NULL,
#'                     pvalues = syn_example_p,
#'                     expression = syn_example_rld,
#'                     p_col_suffix = "pvalue",
#'                     padj_col_suffix = "padj",
#'                     fc_col_suffix = NULL,
#'                     padjust_method = "BH",
#'                     multi_group_prefix = NULL,
#'                     non_sig_name = "Not Significant",
#'                     significance_cutoff = 0.01, 
#'                     fc_cutoff=0.3, 
#'                     label_column = NULL)
#' table(syn_polar@polar$sig) 

polar_coords <- function(sampledata,
                         contrast,
                         pvalues,
                         expression,
                         groups = NULL,
                         p_col_suffix = "pvalues",
                         padj_col_suffix = "padj",
                         fc_col_suffix = NULL,
                         padjust_method = "BH",
                         multi_group_prefix = NULL,
                         non_sig_name = "Not Significant",
                         significance_cutoff = 0.01, 
                         fc_cutoff=0.3, 
                         label_column = NULL, 
                         cutoff_criteria = "pvalue"){
    
    # Check for errors
    if(! cutoff_criteria %in% c("pvalue", "padj")){
      stop("cutoff_criteria must be either 'pvalue', or 'padj' where available")
    }
    if(! class(sampledata) %in% c("data.frame")) {
        stop("sampledata must be a data frame")
    }
    if(! class(expression)[1] %in% c("data.frame", "matrix")) {
        stop("expression must be a data.frame or matrix")
    }
    if(! contrast %in% colnames(sampledata)) {
        stop("contrast is not a column in sampledata")
    }
    if( ! is.null(label_column)){
        if( ! label_column %in% colnames(pvalues)) {
            stop("label_column is not a column in pvalues")
        }
    }
    if(! "ID" %in% colnames(sampledata)) {
        stop("There is no ID column in metadata")
    }
    if(! identical(rownames(expression), rownames(pvalues))){
        stop(paste('The expression row names must be identical to the pvalues', 
                   'row names'))
    }
    if(! identical(colnames(expression), as.character(sampledata$ID))) {
        stop(paste('The expression column names must be identical to the', 
                   'sampledata$ID'))
    }
    if(is.null(label_column)){
        pvalues$label <- rownames(pvalues)
    } else {pvalues$label <- pvalues[, label_column]}
    
    # Ensure groups and contrast column are compatible
    sampledata[, contrast] <- droplevels(factor(sampledata[, contrast]))
    if(length(levels(sampledata[, contrast])) != 3) {
        stop(paste("There number of factors in the comparison column does not",
                   "equal 3"))
    }
    if(any(is.na(sampledata[, contrast]))) {
        stop("There are NAs present in the contrast, please remove these")
    }
    if(is.null(groups)) {groups <- levels(sampledata[, contrast])}
    if(length(groups) != 3) stop("There number of groups does not equal 3")
    if(! all(groups %in% levels(sampledata[, contrast]))) {
        stop('Make sure all groups are in sampledata[, contrast]')
    }
    if(any(grepl("_", groups))) {
        stop("groups must not contain underscores")
    }
    if(any(grepl("_", c(p_col_suffix, fc_col_suffix, padj_col_suffix)))) {
        id <- which(grepl("_", c(p_col_suffix, fc_col_suffix, padj_col_suffix)))
        p <- c("p_col_suffix", "fc_col_suffix", "padj_col_suffix")[id]
        stop(paste(paste(p, collapse=", "), "must not contain underscores"))
    }
    
    comparisons <- c(paste(groups[1], groups[2], sep="_"),
                     paste(groups[2], groups[3], sep="_"),
                     paste(groups[3], groups[1], sep="_"))
    
    # Check column names of correct format exist in pvalues data frame
    found <- num_con <- c()
    invisible(
        lapply(c(p_col_suffix, fc_col_suffix, padj_col_suffix), 
               function(col_suffix){
                   comp_columns <- paste(c(comparisons, multi_group_prefix), 
                                         col_suffix, sep="_")
                   notFinding <- c()
                   if(! all(comp_columns %in% colnames(pvalues))) {
                       notFinding <- comp_columns[!comp_columns %in% 
                                                      colnames(pvalues)]
                       notFinding <- notFinding[!is.na(notFinding)]
                       
                       # check if groups in colnames are reversed
                       check <- strsplit(gsub(paste0("_", col_suffix), "", 
                                              notFinding), "_")
                       
                       for(order_check in check) {
                           og <- paste(order_check[1], order_check[2], 
                                       col_suffix, sep="_")
                           reverse <- paste(order_check[2], order_check[1], 
                                            col_suffix, sep="_")
                           if(reverse %in% colnames(pvalues)){
                               colnames(pvalues)[colnames(pvalues)==
                                                     reverse] <<- og
                               
                               # Need to reverse order for fold change
                               if(! is.null(fc_col_suffix)){
                                   if(col_suffix == fc_col_suffix) {
                                       if(! is.numeric(pvalues[, og])) {
                                           num_con <<- c(num_con, unlist(og))
                                           pvalues[, og] <- 
                                               as.numeric(as.character(
                                                   pvalues[, og]))
                                       }
                                       pvalues[, og] <- -1*pvalues[, og]
                                   }
                               }
                               found <<- c(found, setNames(og, reverse))
                               notFinding <- notFinding[notFinding != og]
                           }
                       }
                   }
                   
                   if(length(notFinding) > 0){  
                       if(length(paste(multi_group_prefix, fc_col_suffix)) > 0 &
                          paste0(multi_group_prefix, "_", fc_col_suffix) %in% 
                          notFinding){
                           notFinding <- 
                               notFinding[notFinding != 
                                              paste0(multi_group_prefix, "_",
                                                     fc_col_suffix)]
                       }
                       if(length(notFinding) > 1) {
                           notFinding <- paste0(
                               "'", paste0(notFinding[1:(length(notFinding)-1)],
                                           collapse="', '"),
                               "' or '", 
                               notFinding[length(notFinding)], "'")
                       }
                       if(length(notFinding) > 0){
                           stop(paste('Cannot find', 
                                      paste0(notFinding, collapse=", "),
                                      'in colnames(pvalues)'))
                       }
                   }
               }))
    
    # Report columns reversed
    if(length(found) > 0){
        found <- unlist(unique(lapply(found, function(x) {
            setNames(paste(unlist(strsplit(x, split="_"))[1:2], collapse="_"),
                     paste(unlist(strsplit(x, split="_"))[2:1], collapse="_"))
        })))
        message(
            paste0("\u2022 ", 
                   paste(sub(",\\s+([^,]+)$", " and \\1", toString(found)), 
                         collapse=", "),
                   ' columns were not found in colnames(pvalues) but ',
                   paste(sub(",\\s+([^,]+)$", " or \\1", 
                             toString(names(found))), collapse=", "),
                   " were. These have been reversed. "))
    }
    
    # extract the statistical parameters
    possible_cols <- paste(
        rep(c(comparisons, multi_group_prefix), 
            each=length(c(p_col_suffix, fc_col_suffix, padj_col_suffix))),
        rep(c(p_col_suffix, fc_col_suffix, padj_col_suffix),
            times = length(c(comparisons, multi_group_prefix))), sep="_")
    possible_cols <- possible_cols[possible_cols %in% colnames(pvalues)]
    
    # If not numeric, convert.
    not_numeric <- names(which(lapply(pvalues[, possible_cols], class) != 
                                   "numeric"))
    if(length(not_numeric) > 0){
    message(paste("\u2022 Some", paste(unique(
        gsub(".*_", "", c(not_numeric, num_con))), collapse=", "), 
        'columns were not originally numeric. These have been converted.'))
    }
    pvalues[, not_numeric] <- vapply(not_numeric, function(x) {
        as.numeric(as.character(pvalues[, x]))
    }, FUN.VALUE=rep(0, nrow(pvalues)))
    if(length(which(colSums(is.na(pvalues)) == nrow(pvalues)))){
        stop(paste(paste(names(which(
            colSums(is.na(pvalues)) == nrow(pvalues))), collapse=", "), 
            'contains only NA or could not be converted to a numeric.'))
    }
    
    pvalues <- pvalues[, c(possible_cols, "label")]
    colnames(pvalues) <- gsub(p_col_suffix, "pvalue", colnames(pvalues))
    
    
    # If adjusted p is not available, calculate
    if(is.null(padj_col_suffix) & ! is.null(padjust_method)) {
        for(comp in c(comparisons, multi_group_prefix)){
            pvalues$new <- p.adjust(pvalues[, paste(comp, p_col_suffix, 
                                                    sep="_")],
                                    method = padjust_method)
            colnames(pvalues)[colnames(pvalues) == "new"] <- 
                paste0(comp, "_padj")
            possible_cols <- c(possible_cols, paste0(comp, "_padj"))
        }
        padj_col_suffix <- "padj"
    } else if(! is.null(padj_col_suffix)){
      colnames(pvalues) <- gsub(padj_col_suffix, "padj", colnames(pvalues))
    }
 
    
    if(! is.null(fc_col_suffix)) {
        colnames(pvalues) <- gsub(fc_col_suffix, "logFC", colnames(pvalues))
    }
    
    if(! is.numeric(significance_cutoff)) {
        stop('significance_cutoff must be a numeric')
    }
    if(! (significance_cutoff >= 0 & significance_cutoff <= 1)) {
        stop('significance_cutoff must be between 0 and 1')
    }
    if(! is.numeric(fc_cutoff)) stop('fc_cutoff must be a numeric')
    
    # Check alignement
    if(identical(as.character(sampledata$ID), colnames(expression)) == FALSE) {
        stop("sampledata and expression data not properly aligned")
    }
    
    # Capture the expression score for each group (fc and z-score)
    expression_scaled <- t(scale(t(expression)))
    polar_colours <- as.data.frame(
        vapply(levels(sampledata[, contrast]), function(x) {
            rowMeans(expression_scaled[, droplevels(sampledata[, contrast])==x])
        }, FUN.VALUE=rep(0, nrow(expression))))
    polar_coloursFC <- as.data.frame(
        vapply(levels(sampledata[, contrast]), function(x) {
            rowMeans(expression[, droplevels(sampledata[, contrast])==x])
        }, FUN.VALUE=rep(0, nrow(expression))))
    
    
    if(! identical(rownames(expression), rownames(pvalues))) {
        stop('expression and pvalues not aligned (rownames are not identical)')
    }
    
    polar_colours$Name <- rownames(polar_colours)
    sampledata$contrast <- droplevels(sampledata[, contrast])
    contrast_groups <- levels(sampledata$contrast)
    
    # Calculate the polar coordinates (uses radians)
    polar_colours$y_zscore <- sinpi(1/3)*(
        polar_colours[, contrast_groups[2]] -
            polar_colours[, contrast_groups[3]])
    polar_colours$x_zscore <- polar_colours[, contrast_groups[1]] -
        (cospi(1/3)*(polar_colours[, contrast_groups[3]] +
                         polar_colours[, contrast_groups[2]]))
    polar_colours$y_fc <- sinpi(1/3)*(polar_coloursFC[,contrast_groups[2]] -
                                          polar_coloursFC[,contrast_groups[3]])
    polar_colours$x_fc <- polar_coloursFC[,contrast_groups[1]] -
        (cospi(1/3)*(polar_coloursFC[,contrast_groups[3]] +
                         polar_coloursFC[,contrast_groups[2]]))
    
    # angle runs from -0.5 to +0.5
    polar_colours$angle <- atan2(polar_colours$y_zscore,
                                 polar_colours$x_zscore)/(2*pi)
    polar_colours$angle_degrees <- (polar_colours$angle %% 1)*360
    
    # rotate and modulus
    polar_colours$angle <- (polar_colours$angle + 2/3) %% 1
    polar_colours$r_zscore <- with(polar_colours, sqrt(x_zscore^2 + y_zscore^2))
    polar_colours$r_fc <- with(polar_colours, sqrt(x_fc^2 + y_fc^2))
    
    # pick the most highly expressed group
    if(! length(which(grepl(cutoff_criteria, colnames(pvalues)))) %in% c(3, 4)){
      stop("cutoff_criteria must have corresponding columns in the pvalues 
           data frame (either pvalue or padj if available) for each comparison")
    }
    groups <- as.character(max.col(polar_colours[, 1:3]))
    if(! is.null(multi_group_prefix)){
        groups[pvalues[, paste(multi_group_prefix, cutoff_criteria, sep="_")] 
               >= significance_cutoff] <- 'grey60'
    }
    polar_colours$max_exp <- colnames(polar_colours)[max.col(
        polar_colours[, 1:3])]
    comp_map <- colnames(polar_colours)[1:3]
    comp_cols <- colnames(pvalues)[grepl(cutoff_criteria, 
                                         colnames(pvalues)) & 
                                       colnames(pvalues) != 
                                       paste(multi_group_prefix, "pvalue")]
    
    groups[pvalues[,comp_cols[1]] >= significance_cutoff &
               pvalues[,comp_cols[2]] >= significance_cutoff &
               pvalues[,comp_cols[3]] >= significance_cutoff] <-
        'grey60'
    polar_colours$max_exp[pvalues[,comp_cols[1]] >= significance_cutoff &
                             pvalues[,comp_cols[2]] >= significance_cutoff &
                             pvalues[,comp_cols[3]] >= significance_cutoff] <-
        non_sig_name
    
    # Calculate which significance group each gene belongs to
    index <- groups != 'grey60' 
    if (any(index != FALSE)){
        pairwise_comp <- data.frame(pvalues[index, comp_cols],
                                    row.names = rownames(pvalues)[index])
        
        # Determine which groups are significant
        pairwise <- data.frame(ifelse(pairwise_comp <= significance_cutoff, 
                                      "1", "0"))
        pairwise$min <- as.numeric(
            apply(polar_colours[index, 1:3], 1, function(x) {
                as.numeric(which.min(x))
            }))
        pairwise$min2 <- comp_map[match(colnames(polar_colours)[1:3],
                                        comp_map)][pairwise$min]
        
        # create a string determiing: min (of ABC),  sig AvB,  sig BvC, sig CvA
        pairwiseSig <- paste0(pairwise$min2, pairwise[,1],
                              pairwise[,2], pairwise[,3])
        sig_res <- c(paste0(comp_map[2], "1..|", comp_map[3], "..1"),
                     paste0(comp_map[1], "1..|", comp_map[3], ".1."),
                     paste0(comp_map[1], "..1|", comp_map[2], ".1."))
        max_res <- c(paste0(comp_map[3], ".11"),
                     paste0(comp_map[1], "1.1"),
                     paste0(comp_map[2], "11."))
        pairwise$max_exp <- non_sig_name
        
        # Fill in colours - for up in one group
        for(iter in 1:3){
            pairwise$max_exp[grep(sig_res[iter], pairwiseSig)] <-
                colnames(polar_colours)[iter]
        }
        
        # Fill colours if up in two groups
        for(iter in 1:3){
            iterIDs <- sort(c(iter, iter%%3 + 1))
            pairwise$max_exp[grep(max_res[iter], pairwiseSig)] <-
                paste(colnames(polar_colours)[iterIDs[1]], "+",
                      colnames(polar_colours)[iterIDs[2]], "+", sep = "")
        }
        
        pairwise$Name <- rownames(pairwise)
        polar_colours$sig <- pairwise$max_exp[match(rownames(polar_colours),
                                                   rownames(pairwise))]
    }
    polar_colours$sig[is.na(polar_colours$sig)] <- non_sig_name
    polar_colours$sig[polar_colours$r_fc < fc_cutoff] <- non_sig_name
    if(! is.null(multi_group_prefix)){
        polar_colours$sig[pvalues[, paste0(multi_group_prefix, "_", 
                                           cutoff_criteria)] >= 
                              significance_cutoff] <- non_sig_name
    } 
    polar_colours$sig <- as.character(polar_colours$sig)
    polar_colours$sig[polar_colours$sig != non_sig_name & 
                          (! grepl('\\+', polar_colours$sig)) ] <- 
        paste0(polar_colours$sig[polar_colours$sig != non_sig_name & 
                                     (! grepl('\\+', polar_colours$sig)) ], '+')
    polar_colours$sig <- factor(polar_colours$sig)
    polar_colours <- polar_colours[, c("Name",
                                       comp_map,
                                       "y_zscore", "x_zscore", "r_zscore",
                                       "x_fc", "y_fc", "r_fc",
                                       "angle",
                                       "angle_degrees",
                                       "max_exp",
                                       "sig")]
    
    if(! identical(rownames(polar_colours), rownames(pvalues))) {
        stop('Misalignment of rows in polar and pvalues')
    }
    
    colnames(polar_colours)[colnames(polar_colours) %in% contrast_groups] <- 
        paste0(colnames(polar_colours)[colnames(polar_colours) %in% 
                                          contrast_groups], "_axis")
    polar_colours$label <- pvalues$label
    
    methods::new("polar",
                 polar = polar_colours,
                 pvalues = pvalues,
                 sampledata = sampledata,
                 contrast = contrast,
                 multi_group_test = multi_group_prefix,
                 expression  = expression,
                 non_sig_name = non_sig_name)
} 




