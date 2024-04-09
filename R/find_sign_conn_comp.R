#' Finding significant connected components
#' @param NRsummary data.frame with NRsummary
#' @param NRpThr threshold for NR p value
#' @param minRank minimum rank
#' @param maxRank maximum rank
#' @return significant connected component
#' @export

find_sign_conn_comp <- function(NRsummary=NULL, NRpThr=0.05, minRank=50, maxRank=300){

    NRsummaryDer <- NRsummary
    #NRsummaryDer$rank <- 1:nrow(NRsummary)
    NRsummaryDer$log10p <- log10(NRsummary$p)
    NRsummaryDer$p_der <- 
        c(NRsummaryDer$log10p[-1] - 
              NRsummaryDer$log10p[-length(NRsummaryDer$log10p)], 0)
    NRsummaryDer$critical <- 0
    NRsummaryDer$selected <- 0
    NRsummaryDer$critical[NRsummaryDer$p <
                           NRpThr & NRsummaryDer$rank >= 
                           minRank & NRsummaryDer$rank <= maxRank] <- 1

    # The critical point with the highest derivative
    critical <- NULL
    if(any(NRsummaryDer$critical==1)){
        selectedRank <- NRsummaryDer[NRsummaryDer$critical==1, ]
        selectedRank <- selectedRank$rank[ which.max(selectedRank$p_der)]
        NRsummaryDer$selected[1:selectedRank] <- 1
        critical <- max(NRsummaryDer$rank[NRsummaryDer$selected==1])
    }

    return(list(signCompTable=NRsummaryDer, critical=critical))

}
