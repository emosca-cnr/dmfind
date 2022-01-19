########################################################################################################
###########################################  INTRODUCTION   ############################################
########################################################################################################
library(dmfind)
library(htmltools)
library(igraph)
library(BiocStyle)
data(X0, W, g400esm)
#######################################################################################################
#####################################      NETWORK DIFFUSION     ######################################
#######################################################################################################
Xs <- ND(X0, W)$Xs
g <- g400esm
jpeg("X0_Xs.jpg", width = 180, height = 240, units="mm", res=300)
plot(X0, Xs, xlab="X0", ylab = "Xs")
dev.off()
#######################################################################################################
#######################################     EPSILON CHOICE     ########################################
#######################################################################################################
eps <- matrix(c(0.01, 1, 10, 100), ncol=1)
out <- dmfind::eval_eps(X0, Xs, eps, g)
dmfind::plot_omega_eps(out)

jpeg("eval_eps_X0.jpg", width = 180, height = 240, units="mm", res=300)
par(mfrow=c(3, 2))
for(i in 1:nrow(eps)){
    plot(X0[, 1], out$S[[i]][, 1], xlab="X0", ylab="S", main = eps[i, ], pch=".")
    top_n <- which(rank(-out$S[[i]][, 1]) <= 300)
    points(X0[top_n, 1], out$S[[i]][top_n, 1], col="red", pch=16, cex=0.7)
}
dev.off() 
#######################################################################################################
##################################     NETWORK SMOOTHED INDEX     #####################################
#######################################################################################################
resS <- calc_adjND(X0 = X0, W = W, eps = 1, k = 49, mode = "S", mc.cores = 2)
scores <- merge(X0, resS$Sp, by=0, sort=F, all.y = TRUE)
rownames(scores) <- scores$Row.names
scores$Row.names <- NULL
names(scores) <- c("X0", "Sp")
init_sig_genes <- scores[order(-scores$X0), ]
dmfind::plot_S(resS, X0, init_sig_genes = init_sig_genes)
#######################################################################################################
####################################     NETWORK RESAMPLING     #######################################
#######################################################################################################
nr_res <- NR(g, sort(resS$Sp[, 1], decreasing = T), k = 99,
             mc.cores = 2)
sigcomp <- find_sign_conn_comp(NR_summary = nr_res$NR_summary, NR_p_thr = 0.2)
plot_NR(nr_res, sign_comp_table = sigcomp$sign_comp_table)
#######################################################################################################
###################################     ENRICHMENT ANALYSIS     #######################################
#######################################################################################################
scores <- scores[order(-scores$X0), ]
compare_X0_ND(g, X0_ranked_by_ND = setNames(scores$X0, rownames(scores)),
              X0_ranked_by_X0 = sort(setNames(scores$X0, rownames(scores)), 
                                     decreasing = T), do.plot = T, file = "cmp_X0_ND.jpg", norm = "n")

ae_res <- assess_enrichment(G=g, top_list = rownames(scores), 
                            ranked_vector_X0 = sort(setNames(scores$X0, rownames(scores)), 
                                                    decreasing = T), do.plot = T, critical = sigcomp$critical,
                            file="top_net_enrichment.jpg")
#######################################################################################################
######################################     TOP NETWORK     ############################################
#######################################################################################################
top_network <- extract_module(graph = g, selected_vertices =
                                  nr_res$NR_summary$id, X0 = X0, 
                              min_subnet_size = 2)
#######################################################################################################
#########################     COMMUNITIES AND TOPOLOGICAL ANALYISIS     ###############################
#######################################################################################################
res_comm <- dmfind::find_communities(top_network)
res_comm$info
res_comm$info$algorithm[which.max(res_comm$info$modularity)]
res_c <- dmfind::assess_centrality(top_network, g)
#######################################################################################################
################    Define nodes roles based on partecipation coefficient and z-score    ##############
#######################################################################################################
V(top_network)$comm_id <- res_comm$comm[[res_comm$info$algorithm[which.max(res_comm$info$modularity)]]]$membership[match(V(top_network)$name, res_comm$comm[[res_comm$info$algorithm[which.max(res_comm$info$modularity)]]]$names)]
V(top_network)$label <- " "
dmfind::pc_wmz(top_network, do_plot=TRUE)

