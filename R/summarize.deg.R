#' Summarize the patient-specific DEGs
#'
#' Summarize the binary differential expression matrix transformed by \code{\link{bi.deg}} or the cross-validated DEGs predicted by \code{\link{deg.specific}}.
#'
#' @param res.deg a list returned by \code{\link{deg.specific}}
#' @param max.n the number of DEGs to display
#' @param ... other setting
#'
#' @author Guofeng Meng
#'
#'
#' @details This function summarizes the gene number before and after cross-validation.
#'
#' @return A data.frame.
#'
#' @examples
#' report.deg <- summarize(deg.spc)
#' @export

summarize.deg.specific <- function(res.deg, max.n = 10, ...) {
    all.pas = length(res.deg[["decd.input"]]$patients)
    cc1 = apply(res.deg[["decd.input"]]$deg, 2, function(x) length(x[x == 1]))
    cc2 = apply(res.deg[["decd.input"]]$deg, 2, function(x) length(x[x == -1]))
    
    pas = names(res.deg)
    pas = pas[pas != "decd.input" & pas != "decd.test"]
    deg.gen.len = vapply(pas, function(x) length(res.deg[[x]]$genes), 10)
    deg.pa.len = vapply(pas, function(x) length(res.deg[[x]]$patients), 10)
    
    pas = names(sort(deg.gen.len, decreasing = TRUE))
    dd1 = vapply(pas, function(x) {
        y = (res.deg[["decd.input"]]$deg)[res.deg[[x]]$genes, x]
        length(y[y == 1])
    }, 1)
    
    dd2 = vapply(pas, function(x) {
        y = (res.deg[["decd.input"]]$deg)[res.deg[[x]]$genes, x]
        length(y[y == -1])
    },1)
    sim = vapply(pas, function(x) res.deg[[x]]$sc, 0.1)
    cutoff = vapply(pas, function(x) res.deg[[x]]$cutoff, 0.1)
    output = data.frame(sample.id = pas, up.regulation.raw = cc1[pas], down.regulation.raw = cc2[pas], 
        neighbor.patients = deg.pa.len[pas], up.regulation.valdiated = dd1[pas], 
        down.regulation.validated = dd2[pas], similarity = sim, cutoff = cutoff)
    output = output[order(output$up.regulation.valdiated + output$down.regulation.validated, 
        decreasing = TRUE), ]
    return(output[seq_len(min(max.n, length(pas))), ])
}

