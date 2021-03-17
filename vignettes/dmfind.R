## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("/path_to_package/dmfind_xxx.tar.gz", repos = NULL)
#  library(dmfind)

## ---- eval=TRUE, message=FALSE------------------------------------------------
library(igraph)
library(dmfind)
data(g400esm, W, X0, package="dmfind")

summary(g400esm)
head(X0)


## ---- eval=TRUE, message=FALSE------------------------------------------------
library(dmfind)
Xs <- ND(X0, W)$Xs
rank_X0 <- rank(-X0, ties.method = "min")
rank_Xs <- rank(-Xs, ties.method = "min")
par(mar=c(4, 4, 2, 1))
par(mfrow=c(1, 2))
L <- layout_nicely(g400esm)
plot(X0[, 1], Xs[, 1], pch=16)
text(X0[, 1], Xs[, 1], V(g400esm)$label[match(rownames(X0), V(g400esm)$name)], pos = 1, cex=0.8)
V(g400esm)$shape <- ifelse(rank(-X0[match(V(g400esm)$name, rownames(Xs)), 1])<=6, "square", "circle")
V(g400esm)$color <- heat.colors(nrow(Xs))[rank_Xs[match(V(g400esm)$name, rownames(X0))]]
par(mar=c(0.1, 0.1, 0.1, 0.1))
plot(g400esm, vertex.size=8, vertex.label.cex=0.5,
     edge.width=1.5, edge.color="hotpink",
     vertex.label=rank_Xs[match(V(g400esm)$name, rownames(X0))],
     vertex.label.cex=1.3,
     layout=L
)


