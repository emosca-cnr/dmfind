#' assess enrichment of the top networks
#' @param G igraph object
#' @param topList ranked list of vertex names that will be used to define top networks
#' @param ranks ranks of topList that will be assessed
#' @param X0Vector named numeric vector that will be tested with GSEA or ORA. In case of GSEA it will be ranked by decreasing orderg, while in the case of ORA the names of the X0 values grater than 0 will be tested for enrichment, while all X0 names will be the universe.
#' @param k number of permutations
#'
#' @export

assess_enrichment <-
  function(G = NULL,
           topList = NULL,
           ranks = NULL,
           X0Vector = NULL,
           type = c("gsea", "ora"),
           k = 99,
           critical = NULL,
           minComponentSize = 2,
           minNetSize = 10,
           minKNes = 10,
           BPPARAMGsl = NULL,
           BPPARAMK = NULL) {
    type <- match.arg(type, c("gsea", "ora"))
    
    if (is.null(ranks)) {
      ranks <- seq(10, length(topList), by = 10)
    }
    en_res <- setNames(vector("list", length(ranks)), ranks)
    net_ccf <- en_res
    
    #defines the increasing top networks to be tested
    for (i in 1:length(ranks)) {
      net_ccf[[i]] <-
        get_nconn_comp(induced_subgraph(G, V(G)$name[V(G)$name %in% topList[1:ranks[i]]]), n =
                         minComponentSize)
    }
    
    gsl <- lapply(net_ccf, function(i_net)
      V(i_net)$name)
    gsl <- gsl[lengths(gsl) >= minNetSize]
    net_ccf <- net_ccf[names(net_ccf) %in% names(gsl)]
    
    if (length(gsl) > 0) {
      #check X0 to choose ORA or GSEA
      if (type == "ora") {
        cat("ORA\n")
        
        wb <- names(X0Vector)[X0Vector > 0]
        if (length(wb) > 0) {
          bb <- names(X0Vector)[!names(X0Vector) %in% wb]
          en_res <- ora(
            wb = wb,
            bb = bb,
            gsl = gsl,
            p_adj_method = 'fdr'
          )
        } else{
          stop("Can' find any element of X0Vector > 0\n.")
        }
        
      } else{
        cat("GSEA\n")
        
        #if there are zeros... a noise is added to avoid ties in the ranked gene list
        if (sum(X0Vector == 0) > 1) {
          cat("Detected multiple null values in X0. Adding a noise to avoid ties.\n")
          idx_zeros <- which(X0Vector == 0)
          min_pos <- min(X0Vector[-idx_zeros])
          X0Vector[idx_zeros] <-
            abs(rnorm(
              n = length(idx_zeros),
              mean = 0,
              sd = min_pos * 10 ^ -6
            ))
          print(summary(X0Vector[-idx_zeros]))
          print(summary(X0Vector[idx_zeros]))
        }
        
        en_res <-
          gsea(
            rl = matrix(
              X0Vector,
              ncol = 1,
              dimnames = list(names(X0Vector), "X0")
            ),
            gsl = gsl,
            k = k,
            minKNes = minKNes,
            BPPARAMGsl = BPPARAMGsl,
            BPPARAMK = BPPARAMK
          )
        en_res <- en_res$gs_table$X0
        
      }
      
      #en_res_df <- do.call(rbind, lapply(en_res, function(x) x$gs_table$X0))
      #en_res_df$id <- ranks
      en_res$size <- unlist(lapply(net_ccf, function(x)
        vcount(x)))
      
      colnames(en_res) <- replace(colnames(en_res), colnames(en_res) == "id", "rank")
      
      return(list(en_summary = en_res, net_ccf = net_ccf))
      
    } else{
      cat("No networks with at least ", minNetSize, ".\n")
      return(NULL)
    }
  }
