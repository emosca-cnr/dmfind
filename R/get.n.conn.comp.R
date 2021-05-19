get.n.conn.comp <- function(x, n){
  conn <- clusters(x)
  clstr.ok <- which(conn$csize >= n)
  nodes.ok <- which(conn$membership %in% clstr.ok)
  return(induced.subgraph(x, nodes.ok))
}
