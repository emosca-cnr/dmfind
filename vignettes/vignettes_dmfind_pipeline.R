setwd("/home/bioinformatics/vnale/dmfind/vignettes/imgs_2")
########################################################################################################
###########################################  INTRODUCTION   ############################################
########################################################################################################
library(dmfind)
library(igraph)
data(X0, W, g400esm)
#######################################################################################################
#####################################      NETWORK DIFFUSION     ######################################
#######################################################################################################
Xs <- ND(X0, W)$Xs
g <- g400esm
jpeg("X0_Xs.jpg", width = 180, height = 180, units="mm", res=300)
plot(X0, Xs, xlab="X0", ylab = "Xs")
dev.off()
#######################################################################################################
#######################################     EPSILON CHOICE     ########################################
#######################################################################################################
eps <- matrix(c(0.01, 1, 10, 100), ncol=1)
out <- dmfind::eval_eps(X0, Xs, eps, g)
dmfind::plot_omega_eps(out)
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
                                                    decreasing = T), do.plot = T, critical = 40,
                            file="top_net_enrichment.jpg")
#######################################################################################################
######################################     TOP NETWORK     ############################################
#######################################################################################################
top_network <- extract_module(graph = g, selected_vertices =
                                  nr_res$NR_summary$id[1:40], X0 = X0,
                              min_subnet_size = 2)
#######################################################################################################
#########################     COMMUNITIES AND TOPOLOGICAL ANALYISIS     ###############################
#######################################################################################################
res_comm <- dmfind::find_communities(top_network)
res_comm$info
# algorithm modularity n
# 1 fastgreedy  0.4399270 6
# 2    labprop  0.1974721 4
# 3   walktrap  0.4445886 7
# 4      eigen  0.4264225 6
# 5   multilev  0.4562188 5 <-
# 6    infomap  0.4538158 6
res_comm$info$algorithm[which.max(res_comm$info$modularity)]
res_c <- dmfind::assess_centrality(top_network, g)

pal <- rainbow(length(unique(V(top_network)$comm_id)))
lo <- plot_network(top_network, color_by = "comm_id", label_by = "", pal = pal, vertex.size=6,
                   color_quant=F, community = res_comm$comm[[res_comm$info$algorithm[which.max(res_comm$info$modularity)]]],
                   comm_w_in=2, comm_w_b = 1, plot_outfile = "comm_network.png")

#######################################################################################################
################    Define nodes roles based on partecipation coefficient and z-score    ##############
#######################################################################################################
V(top_network)$comm_id <- res_comm$comm[[res_comm$info$algorithm[which.max(res_comm$info$modularity)]]]$membership[match(V(top_network)$name, res_comm$comm[[res_comm$info$algorithm[which.max(res_comm$info$modularity)]]]$names)]
V(top_network)$label <- " "
sym <- as.data.frame(org.Hs.egSYMBOL)
V(top_network)$label <- sym$symbol[match(V(top_network)$name, sym$gene_id)]

dmfind::pc_wmz(top_network, do_plot=TRUE, all_names = T)

#######################################################################################################
#####################################    Communities network    #######################################
#######################################################################################################
commNet <- comm_net(top_network)
plot_network(commNet[[1]], color_by = "comm_id", pal=pal, comm_w_in = 0.1, comm_w_b = 10,
             plot_outfile = "comm_net.jpg", legend.off = T, min_subnet_size = 2)

#######################################################################################################
######################################    Network Comparison   ########################################
#######################################################################################################
n1 <- as.character(sample(1:70, 45))
n2 <- as.character(sample(1:70, 52))
n3 <- as.character(sample(1:70, 35))

g1 <- barabasi.game(45, directed = F)
V(g1)$name <- n1
l1 <- paste0("sym", n1)
V(g1)$label <- l1

g2 <- barabasi.game(52, directed = F)
V(g2)$name <- n2
l2 <- paste0("sym", n2)
V(g2)$label <- l2

g3 <- barabasi.game(35, directed = F)
V(g3)$name <- n3
l3 <- paste0("sym", n3)
V(g3)$label <- l3


g <- list(g1, g2, g3)
netcom <- NetComp(g, union=T, intersection=T, centrMeas = T, ji = T, oc = T)
