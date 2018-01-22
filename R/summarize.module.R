#' Summarize the DEG modules
#'
#' @param res.module a 'seed.module' object returned by \code{\link{seed.module}}
#' @param max.n the number of modules to display
#' @param ... other setting
#'
#' @author Guofeng Meng
#' @references
#'
#' @details This function summarize the DEG modules
#' @return A data.frame.
#'
#' @examples
#' report.module <- summarize(seed.mod)
#' @export

summarize.seed.module <- function(res.module, max.n = 10, ...) {
    all.pas = colnames(res.module[["decd.input"]]$deg)
    aa = apply(res.module[["decd.input"]]$deg, 2, function(x) length(x[x != 0]))
    pas = names(res.module)
    pas = pas[pas != "decd.input" & pas != "decd.specific" & pas != "decd.clustering"]
    deg.gen.len1 = vapply(pas, function(x) length(res.module[[x]][["max.genes"]][["genes"]]), 10)
    deg.pa.len1 = vapply(pas, function(x) length(res.module[[x]][["max.genes"]][["patients"]]), 10)
    sc1 = unlist(vapply(pas, function(x) res.module[[x]][["max.genes"]][["sc"]], 0.1))
    deg.gen.len2 = vapply(pas, function(x) length(res.module[[x]][["model"]][["genes"]]), 10)
    deg.pa.len2 = vapply(pas, function(x) length(res.module[[x]][["model"]][["patients"]]), 10)
    sc2 = unlist(vapply(pas, function(x) res.module[[x]][["model"]][["sc"]],0.1))
    deg.gen.len3 = vapply(pas, function(x) length(res.module[[x]][["max.patients"]][["genes"]]), 10)
    deg.pa.len3 = vapply(pas, function(x) length(res.module[[x]][["max.patients"]][["patients"]]), 10)
    sc3 = unlist(vapply(pas, function(x) res.module[[x]][["max.patients"]][["sc"]], 0.1))
    overlap = vapply(pas, function(x) length(res.module[[x]][["max.patients"]][["cutoff"]]),0.9)
    
    output = data.frame(modules = pas, samples = rep(length(all.pas), length(pas)), 
        genes.max.genes = deg.gen.len1, patients.max.genes = deg.pa.len1, sc.max.genes = sc1, 
        genes.model = deg.gen.len2, patients.model = deg.pa.len2, sc.model = sc2, 
        genes.max.patients = deg.gen.len3, patients.max.patients = deg.pa.len3, sc.max.patients = sc3, 
        overlap = overlap)
    # names(output)<-c('modules', 'samples', 'genes(max.genes)',
    # 'patients(max.genes)', 'sc(max.genes)', 'genes(model)','patients(model)',
    # 'sc(model)', 'genes(max.patients)', 'patients(max.patients)',
    # 'sc(max.patients)', 'overlap')
    return(output[1:min(max.n, dim(output)[1])])
}


#' @rdname summarize.seed.module
#' @export
summarize <- function(...) {
    UseMethod("summarize")
}


#' Summarize the DEG modules
#'
#' @name summarize.cluster.module
#'
#' @param res.module a 'cluster.module' object returned by \code{\link{cluster.module}}
#' @param max.n the number of modules to display
#' @param ... other setting
#'
#' @author Guofeng Meng
#' @references
#'
#' @details This function summarize the DEG modules
#' @return A data.frame.
#'
#' @examples
#' 
#' report.module <- summarize(cluster.mod)
#' 
#' @export
summarize.cluster.module <- function(res.module, max.n = 10, ...) {
    all.pas = colnames(res.module[["decd.input"]]$deg)
    aa = apply(res.module[["decd.input"]]$deg, 2, function(x) length(x[x != 0]))
    pas = names(res.module)
    pas = pas[pas != "decd.input" & pas != "decd.speicific" & pas != "decd.clustering"]
    deg.gen.len1 = vapply(pas, function(x) length(res.module[[x]][["max.genes"]][["genes"]]),10)
    deg.pa.len1 = vapply(pas, function(x) length(res.module[[x]][["max.genes"]][["patients"]]),10)
    sc1 = unlist(vapply(pas, function(x) res.module[[x]][["max.genes"]][["sc"]], 0.9))
    deg.gen.len2 = vapply(pas, function(x) length(res.module[[x]][["model"]][["genes"]]),10)
    deg.pa.len2 = vapply(pas, function(x) length(res.module[[x]][["model"]][["patients"]]),10)
    sc2 = unlist(vapply(pas, function(x) res.module[[x]][["model"]][["sc"]],0.9))
    deg.gen.len3 = vapply(pas, function(x) length(res.module[[x]][["max.patients"]][["genes"]]),10)
    deg.pa.len3 = vapply(pas, function(x) length(res.module[[x]][["max.patients"]][["patients"]]),10)
    sc3 = unlist(vapply(pas, function(x) res.module[[x]][["max.patients"]][["sc"]],0.1))
    overlap = vapply(pas, function(x) length(res.module[[x]][["max.patients"]][["cutoff"]]),10)
    
    output = data.frame(modules = pas, samples = length(all.pas), genes.max.genes = deg.gen.len1, 
        patients.max.genes = deg.pa.len1, sc.max.genes = sc1, genes.model = deg.gen.len2, 
        patients.model = deg.pa.len2, sc.model = sc2, genes.max.patients = deg.gen.len3, 
        patients.max.patients = deg.pa.len3, sc.max.patients = sc3, overlap = overlap)
    
    # names(output)<-c('modules', 'samples', 'genes(max.genes)',
    # 'patients(max.genes)', 'sc(max.genes)', 'genes(model)','patients(model)',
    # 'sc(model)', 'genes(max.patients)', 'patients(max.patients)',
    # 'sc(max.patients)', 'overlap')
    return(output[1:min(max.n, dim(output)[1]), ])
}

