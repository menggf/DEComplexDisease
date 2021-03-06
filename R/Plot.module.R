#' Plot the DEGs modules shared by patients
#'
#' Plot the DEGs modules
#'
#' @import ComplexHeatmap
#' @import grid
#'
#' @param res.module a 'seed.module' object returned by \code{\link{seed.module}}
#' @param ann a data.frame for the patient annotation
#' @param deg a 'deg' to display. It is returned by \code{\link{bi.deg}}
#' @param col.order the order of column in heatmap
#' @param show.mods a vector, the modules to display
#' @param overlap the similarity cutoff to display as carrying the module
#' @param dissimilarity the similarity cutoff to display as not carrying the module
#' @param max.n the maximum number of modules to display
#' @param type the module type to display
#' @param label.col the color to label
#' @param ... other setting of 'oncoPrint'
#'
#' @author Guofeng Meng
#'
#'
#' @details This function is to display the relationship of the predicted DEG modules and the patients.
#'
#' 'deg' can be set to display the modules from different datasets, e.g. the modules predicted from disease A and display them in the binary DEG matrix of disease B.
#'
#' The output is a heatmap Plot where the modules with maximum observations are showed.
#' @return A heatmap plot
#'
#' @examples
#' Plot(seed.mod, ann.er, max.n=5)
#' Plot(seed.mod, ann.er, deg=deg, max.n=5)
#' @export

Plot.seed.module <- function(res.module, ann = NULL, deg = NULL, col.order = NULL, 
    show.mods = NULL, overlap = NULL, dissimilarity = NULL, max.n = min(length(res.module), 
        30), type = c("model", "max.patients", "max.genes")[1], label.col = "blue", 
    ...) {
    mods = names(res.module)
    mods = mods[mods != "decd.specific" & mods != "decd.input" & mods != "decd.clustering"]
    if (!any(c("model", "max.patients", "max.genes") == type)) 
        stop("Error: type: should one of model, max.patients and max.genes!")
    if (length(mods) == 0) 
        stop("Error: No module is available")
    if (is.null(show.mods)) {
        show.mods = .select.mod(res.module, max.n, type = type)
    } else {
        show.mods = show.mods[show.mods %in% mods]
    }
    if (length(show.mods) <= 1) 
        stop("Error: show.mods: no id is recognized")
    
    if (is.null(deg)) {
        deg = res.module[["decd.input"]][["deg"]]
    }
    
    
    pas = colnames(deg)
    ges = row.names(deg)
    if (is.null(overlap)) 
        overlap = res.module[["decd.input"]][["overlap"]]
    
    if (is.null(dissimilarity)) {
        mycols = label.col
        ck = overlap
        mylabs = paste(">= ", ck, sep = "")
    } else {
        mycols = c(label.col, "yellow")
        ck = c(overlap, dissimilarity)
        mylabs = paste(c(">= ", "< "), ck, sep = "")
    }
    names(mycols) <- mylabs
    
    alter_fun = list()
    alter_fun[["background"]] = function(x, y, w, h) {
        grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = "grey", 
            col = NA))
    }
    alter_fun[[mylabs[1]]] = function(x, y, w, h, col = mycols) {
        grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = label.col, 
            col = NA))
    }
    if (!is.null(dissimilarity)) 
        alter_fun[[mylabs[2]]] = function(x, y, w, h, col = mycols) {
            grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = "yellow", 
                col = NA))
        }
    
    
    mat.shared = t(vapply(show.mods, function(x) {
        seed = res.module[[x]][["seed"]]
        if (x == "M0") {
            # used.ges=res.module[[x]][['max.patients']][['genes']] #module genes
            used.ges = res.module[[x]][[type]][["genes"]]
        } else {
            used.ges = res.module[[x]][[type]][["genes"]]
        }
        
        if (length(used.ges[used.ges %in% ges])/length(used.ges) < 0.5) {
            stop("Error: 'deg': >50% of module genes are not observed")
        }
        if (length(used.ges[used.ges %in% ges])/length(used.ges) < 0.8) {
            warnings("'deg': >20% of module genes are not observed")
        }
        used.seed = seed[used.ges]
        sub.deg = deg[used.ges, ]
        sims = apply(sub.deg, 2, function(z) length(which(z == used.seed))/length(used.seed))
        rr = rep("", length(pas))
        rr[sims >= overlap] = mylabs[1]
        if (!is.null(dissimilarity)) 
            rr[sims < dissimilarity] = mylabs[2]
        
        return(rr)
    }, rep("", length(pas))))
    len = vapply(show.mods, function(x) {
        if (is.null(res.module[[x]][[type]])) 
            return(0)
        if (x == "M0") {
            # return(length(res.module[[x]][['max.patients']][['genes']]));
            return(length(res.module[[x]][[type]][["genes"]]))
        } else {
            return(length(res.module[[x]][[type]][["genes"]]))
        }
    }, 1)
    ha = NULL
    if (!is.null(ann)) {
        has.pas = row.names(ann)
        if (length(which(has.pas %in% pas)) < 0.6 * length(pas)) 
            warnings("Warning: ann: Too few patients has annotation")
        if (length(which(has.pas %in% pas)) < 0.3 * length(pas)) 
            stop("Error: ann: Too few patients has annotation")
        ann.all = as.vector(as.matrix(ann))
        ann.all = ann.all[!is.na(ann.all)]
        cl = rainbow(length(unique(ann.all)))
        names(cl) <- unique(ann.all)
        col.list = lapply(names(ann), function(x) {
            has.ann = unique(as.vector(as.matrix(ann[, x])))
            has.ann = has.ann[!is.na(has.ann)]
            return(cl[has.ann])
        })
        names(col.list) <- names(ann)
        if (dim(ann)[2] == 1) {
            new.ann = as.data.frame(ann[pas, ])
            row.names(new.ann) <- pas
            names(new.ann) <- names(ann)
        } else {
            new.ann = ann[pas, ]
        }
        ha = HeatmapAnnotation(df = new.ann, annotation_height = 0.2, name = names(ann), 
            col = col.list)
    }
    row.names(mat.shared) <- paste(show.mods, "(", len, ")", sep = "")
    colnames(mat.shared) <- pas
    
    if (!is.null(ha)) {
        if (is.null(col.order)) {
            if (is.null(dissimilarity)) {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  bottom_annotation = ha, alter_fun = alter_fun, col = mycols, column_title = "", 
                  show_heatmap_legend = FALSE)
            } else {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  bottom_annotation = ha, alter_fun = alter_fun, col = mycols, column_title = "", 
                  show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Overlap"), 
                  show_pct = FALSE)
            }
        } else {
            if (is.null(dissimilarity)) {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  column_order = col.order, bottom_annotation = ha, alter_fun = alter_fun, 
                  col = mycols, column_title = "", show_heatmap_legend = FALSE, ...)
            } else {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  column_order = col.order, bottom_annotation = ha, alter_fun = alter_fun, 
                  col = mycols, column_title = "", show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Overlap"), 
                  show_pct = FALSE, ...)
            }
        }
    } else {
        if (is.null(col.order)) {
            if (is.null(dissimilarity)) {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  alter_fun = alter_fun, col = mycols, column_title = "", show_heatmap_legend = FALSE, 
                  ...)
            } else {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  alter_fun = alter_fun, col = mycols, column_title = "", show_heatmap_legend = TRUE, 
                  heatmap_legend_param = list(title = "Overlap"), show_pct = FALSE, 
                  ...)
            }
        } else {
            if (is.null(dissimilarity)) {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  column_order = col.order, alter_fun = alter_fun, col = mycols, 
                  column_title = "", show_heatmap_legend = FALSE, ...)
            } else {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  column_order = col.order, alter_fun = alter_fun, col = mycols, 
                  column_title = "", show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Overlap"), 
                  show_pct = FALSE, ...)
            }
        }
    }
}


#' Plot the DEGs modules
#'
#' Plot the DEGs modules
#'
#' @import ComplexHeatmap
#' @import grid
#' @param res.module a 'cluster.module' object returned by  \code{\link{cluster.module}}
#' @param ann a data.frame for the patient annotation
#' @param deg a 'deg' to display. It is returned by \code{\link{bi.deg}}
#' @param col.order the order of column in heatmap
#' @param show.mods a vector, the modules to display
#' @param overlap the similarity cutoff to display as carrying the module
#' @param dissimilarity the similarity cutoff to display as not carrying the module
#' @param max.n the maximum number of modules to display
#' @param type the module type to display
#' @param label.col the color to label
#' @param ... other setting of 'oncoPrint'
#'
#' @author Guofeng Meng
#'
#' @references
#'
#' @details This function is to dispaly the relationship of the predicted DEG modules and the patients.
#'
#' 'deg' can be set to display the modules from different datasets, e.g. the modules predicted from disease A and display them in the binary DEG matrix of disease B.
#'
#' The output is a heatmap Plot where the modules with maximum observations are showed.
#' @return A heatmap plot
#'
#' @examples
#' Plot(cluster.mod, ann.er, max.n=5)
#' Plot(cluster.mod, ann.er, deg=deg, max.n=5)
#' @export
#'
Plot.cluster.module <- function(res.module, ann = NULL, deg = NULL, col.order = NULL, 
    show.mods = NULL, overlap = NULL, dissimilarity = NULL, max.n = min(length(res.module), 
        30), type = c("model", "max.patients", "max.genes")[1], label.col = "blue", 
    ...) {
    mods = names(res.module)
    mods = mods[mods != "decd.specific" & mods != "decd.input" & mods != "decd.clustering"]
    if (!any(c("model", "max.patients", "max.genes") == type)) 
        stop("Error: type: should one of model, max.patients and max.genes!")
    if (length(mods) == 0) 
        stop("Error: No module is available")
    if (is.null(show.mods)) {
        show.mods = .select.mod(res.module, max.n, type = type)
    } else {
        show.mods = show.mods[show.mods %in% mods]
    }
    if (length(show.mods) <= 1) 
        stop("Error: show.mods: no id is recognized")
    
    if (is.null(deg)) {
        deg = res.module[["decd.input"]][["deg"]]
        if (!is.null(res.module[["decd.input"]][["test.patients"]])) 
            deg = deg[, res.module[["decd.input"]][["test.patients"]]]
    }
    
    pas = colnames(deg)
    ges = row.names(deg)
    if (is.null(overlap)) 
        overlap = res.module[["decd.input"]][["overlap"]]
    
    if (is.null(dissimilarity)) {
        mycols = label.col
        ck = overlap
        mylabs = paste(">= ", ck, sep = "")
    } else {
        mycols = c(label.col, "yellow")
        ck = c(overlap, dissimilarity)
        mylabs = paste(c(">= ", "< "), ck, sep = "")
    }
    names(mycols) <- mylabs
    
    alter_fun = list()
    alter_fun[["background"]] = function(x, y, w, h) {
        grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = "grey", 
            col = NA))
    }
    alter_fun[[mylabs[1]]] = function(x, y, w, h, col = mycols) {
        grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = label.col, 
            col = NA))
    }
    if (!is.null(dissimilarity)) 
        alter_fun[[mylabs[2]]] = function(x, y, w, h, col = mycols) {
            grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = "yellow", 
                col = NA))
        }
    
    
    mat.shared = t(vapply(show.mods, function(x) {
        seed = res.module[[x]][["seed"]]
        if (x == "M0") {
            # used.ges=res.module[[x]][['max.patients']][['genes']] #module genes
            used.ges = res.module[[x]][[type]][["genes"]]
        } else {
            used.ges = res.module[[x]][[type]][["genes"]]
        }
        
        if (length(used.ges[used.ges %in% ges])/length(used.ges) < 0.5) {
            stop("Error: 'deg': >50% of module genes are not observed")
        }
        if (length(used.ges[used.ges %in% ges])/length(used.ges) < 0.8) {
            warnings("'deg': >20% of module genes are not observed")
        }
        used.seed = seed[used.ges]
        sub.deg = deg[used.ges, ]
        # sims=apply(sub.deg, 2 , function(z) length(which( z ==
        # used.seed)))/length(used.seed);
        sims = colSums(sub.deg == used.seed)/length(used.seed)
        rr = rep("", length(pas))
        rr[sims >= overlap] = mylabs[1]
        if (!is.null(dissimilarity)) 
            rr[sims < dissimilarity] = mylabs[2]
        
        return(rr)
    },rep("", length(pas))))
    len = vapply(show.mods, function(x) {
        if (is.null(res.module[[x]][[type]])) 
            return(0)
        if (x == "M0") {
            # return(length(res.module[[x]][['max.patients']][['genes']]));
            return(length(res.module[[x]][[type]][["genes"]]))
        } else {
            return(length(res.module[[x]][[type]][["genes"]]))
        }
    }, 1)
    ha = NULL
    
    if (!is.null(ann)) {
        has.pas = row.names(ann)
        if (length(which(has.pas %in% pas)) < 0.6 * length(pas)) 
            warnings("Warning: ann: Too few patients has annotation")
        if (length(which(has.pas %in% pas)) < 0.3 * length(pas)) 
            stop("Error: ann: Too few patients has annotation")
        ann.all = as.vector(as.matrix(ann))
        ann.all = ann.all[!is.na(ann.all)]
        cl = rainbow(length(unique(ann.all)))
        names(cl) <- unique(ann.all)
        col.list = lapply(names(ann), function(x) {
            has.ann = unique(as.vector(as.matrix(ann[, x])))
            has.ann = has.ann[!is.na(has.ann)]
            return(cl[has.ann])
        })
        names(col.list) <- names(ann)
        if (dim(ann)[2] == 1) {
            new.ann = as.data.frame(ann[pas, ])
            row.names(new.ann) <- pas
            names(new.ann) <- names(ann)
        } else {
            new.ann = ann[pas, ]
        }
        ha = HeatmapAnnotation(df = new.ann, annotation_height = 0.2, name = names(ann), 
            col = col.list)
    }
    row.names(mat.shared) <- paste(show.mods, "(", len, ")", sep = "")
    colnames(mat.shared) <- pas
    
    if (!is.null(ha)) {
        if (is.null(col.order)) {
            if (is.null(dissimilarity)) {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  bottom_annotation = ha, alter_fun = alter_fun, col = mycols, column_title = "", 
                  show_heatmap_legend = FALSE)
            } else {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  bottom_annotation = ha, alter_fun = alter_fun, col = mycols, column_title = "", 
                  show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Overlap"), 
                  show_pct = FALSE)
            }
        } else {
            if (is.null(dissimilarity)) {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  column_order = col.order, bottom_annotation = ha, alter_fun = alter_fun, 
                  col = mycols, column_title = "", show_heatmap_legend = FALSE, ...)
            } else {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  column_order = col.order, bottom_annotation = ha, alter_fun = alter_fun, 
                  col = mycols, column_title = "", show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Overlap"), 
                  show_pct = FALSE, ...)
            }
        }
    } else {
        if (is.null(col.order)) {
            if (is.null(dissimilarity)) {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  alter_fun = alter_fun, col = mycols, column_title = "", show_heatmap_legend = FALSE, 
                  ...)
            } else {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  alter_fun = alter_fun, col = mycols, column_title = "", show_heatmap_legend = TRUE, 
                  heatmap_legend_param = list(title = "Overlap"), show_pct = FALSE, 
                  ...)
            }
        } else {
            if (is.null(dissimilarity)) {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  column_order = col.order, alter_fun = alter_fun, col = mycols, 
                  column_title = "", show_heatmap_legend = FALSE, ...)
            } else {
                oncoPrint(mat.shared, get_type = function(x) strsplit(x, ";")[[1]], 
                  column_order = col.order, alter_fun = alter_fun, col = mycols, 
                  column_title = "", show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Overlap"), 
                  show_pct = FALSE, ...)
            }
        }
    }
}
