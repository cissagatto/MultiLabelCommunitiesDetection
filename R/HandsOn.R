# workspace directory
setwd("~/Hands-On-Igraph-R")

library(igraph)
library(tidyr)
library(plyr)
library(dplyr)

##############################################################################
# GETTING THE LABELS
##############################################################################
names_labes = read.csv("~/Hands-On-Igraph-R/data/birds-dataset/NamesLabels/birds-NamesLabels.csv")


##############################################################################
# GETTING THE SIMILARITIES AND WEIGHTS                                       #
##############################################################################
df_jaccard = read.csv("~/Hands-On-Igraph-R/data/birds-similarities/Split-1/jaccard-2.csv")
rownames(df_jaccard) = df_jaccard$X
df_jaccard = df_jaccard[,-1]

df_cd = read.csv("~/Hands-On-Igraph-R/data/birds-similarities/Split-1/birds-conditional-probabilities.csv")
rownames(df_cd) = df_cd$X
df_cd = df_cd[,-1]


##############################################################################
# RÓTULOS ADICIONAIS
##############################################################################
create.pivot.labels <- function(names_labes){
  i = 1
  j = 1
  a = 1
  labels.pivot = c()
  for(i in 1:19){
    for(j in 1:19){
      labels.pivot[a] = names_labes$Labels[i]
      a = a + 1
    }
  }
  return(labels.pivot)
}
from = create.pivot.labels(names_labes)


##############################################################################
# BUILD THE DATA FRAME                                                       #
##############################################################################
df_graph = data.frame(pivot_longer(df_jaccard, cols = 1:19))
colnames(df_graph) = c("to", "jaccard")

weigths_total = data.frame(pivot_longer(df_cd, cols = 1:19))
names(weigths_total)[2] = "weights"

df_graph = cbind(from, df_graph, weigths = weigths_total$weights)


##############################################################################
# STRATIFICATION WITH TRESHOLD                                               #
##############################################################################
df_graph_tr = arrange(df_graph, desc(jaccard))
nrow(df_graph_tr)
df_graph_2 = df_graph_tr

df_graph_tr = df_graph_tr[c(-1:-19),]
nrow(df_graph_tr)

df_graph_tr = df_graph_tr[!(df_graph_tr$jaccard == 0 | df_graph_tr$jaccard < 0.01),]
nrow(df_graph_tr)


##############################################################################
# STRATIFICATION WITH KNN                                                    #
##############################################################################

df_graph_2 = arrange(df_graph_2, desc(from))
nrow(df_graph_2)

df_graph_knn = data.frame()
i = 1
for(i in 1:19){
  # filtra pelo rótulo específico
  res1 = df_graph_2[df_graph_2$from == toString(names_labes$Labels[i]),]
  
  # ordena os pesos do maior para o menor
  knn = res1[order(res1[,4], decreasing=TRUE), ]
  
  # faz jaccard = 0 a partir da linha 5
  knn[5:19,]$jaccard = 0
  
  # monta o data frame final
  df_graph_knn = rbind(df_graph_knn, knn)
  gc()
}


##############################################################################
# BUILD THE GRAPH                                                            #
##############################################################################
jaccard_graph_tr = graph_from_data_frame(df_graph_tr, directed=F) 
jaccard_graph_tr
V(jaccard_graph_tr)
E(jaccard_graph_tr)

jaccard_graph_knn = graph_from_data_frame(df_graph_knn, directed=F) 
jaccard_graph_knn
V(jaccard_graph_knn)
E(jaccard_graph_knn)

par(mfrow=c(1,2))
plot(jaccard_graph_tr)
title(cex.main= 0.8, main = "Treshold Graph")
plot(jaccard_graph_knn)
title(cex.main= 0.8, main = "Knn Graph")


##############################################################################
# COMMUNITIES DETECTION METHODS                                              #
##############################################################################


##############################################################################
# 1. Edge Betweenness (aresta intermediária)
##############################################################################

eb_tr = cluster_edge_betweenness(jaccard_graph_tr, 
                                 weights = jaccard_graph_tr$weights,
                                 edge.betweenness = TRUE,
                                 merges = TRUE, 
                                 bridges = TRUE,
                                 modularity = TRUE, 
                                 membership = TRUE)


eb_knn = cluster_edge_betweenness(jaccard_graph_knn, 
                                  weights = jaccard_graph_tr$weights,
                                  edge.betweenness = TRUE,
                                  merges = TRUE, 
                                  bridges = TRUE,
                                  modularity = TRUE, 
                                  membership = TRUE)

plotar <- function(comunidade, grafo, titulo){
  plot(comunidade, grafo, vertex.size=18, edge.arrow.size=.5, 
       vertex.color="gold", vertex.size=12, 
       vertex.shape = "sphere", vertex.frame.color="orange", 
       vertex.label.color="black", vertex.label.cex=0.7, 
       vertex.label.family = "Times", vertex.label.font = 2, 
       edge.color = "gray", edge.width = 0.5) 
  title(cex.main= 0.8, main = titulo)
}

par(mfrow=c(1,2))
plotar(eb_tr, jaccard_graph_tr, "Edge Betweenness threshold")
plotar(eb_knn, jaccard_graph_knn, "Edge Betweenness knn")

eb_tr
eb_knn

compare(eb_tr, eb_knn)
compare(membership(eb_tr), membership(eb_knn))

eb_tr$removed.edges
eb_knn$removed.edges

eb_tr$edge.betweenness
eb_knn$edge.betweenness

eb_tr$merges
eb_knn$merges

eb_tr$bridges
eb_knn$bridges

eb_tr$membership
eb_knn$membership 
membership(eb_tr)
membership(eb_knn)

sizes(eb_tr)
sizes(eb_knn)

eb_tr$vcount
eb_knn$vcount 

eb_tr$modularity
modularity(eb_tr)
modularity_matrix(jaccard_graph_tr)

eb_knn$modularity
modularity(eb_knn)
modularity_matrix(jaccard_graph_knn)

communities(eb_tr)
communities(eb_knn) 

is_hierarchical(eb_tr)
dend_eb = as.dendrogram(eb_tr)
hc_eb = as.hclust(eb_tr)
hc_eb$order

is_hierarchical(eb_knn)
dend_eb_knn = as.dendrogram(eb_knn)
hc_eb_knn = as.hclust(eb_knn)
hc_eb_knn$order

par(mfrow=c(2,1))
plot_dendrogram(eb_tr)
title(cex.main= 0.8, main = "Dendrogram threshold")
plot_dendrogram(eb_knn)
title(cex.main= 0.8, main = "Dendrogram knn")


##############################################################################
# 2. SPINGLASS
##############################################################################

sc_tr = spinglass.community(jaccard_graph_tr, 
                            weights = jaccard_graph_tr$weights,
                            implementation = "neg", 
                            update.rule = "simple")

sc_knn = spinglass.community(jaccard_graph_knn,
                             weights = jaccard_graph_knn$weights,
                             implementation = "neg", 
                             update.rule = "simple")

par(mfrow=c(1,2))
plotar(sc_tr, jaccard_graph_tr, "Spinglass threshold")
plotar(sc_knn, jaccard_graph_knn, "Spinglass knn")

sc_tr
membership(sc)
sizes(sc_tr)
modularity(sc_tr)
communities(sc_tr)
is_hierarchical(sc_tr)


##############################################################################
# 3. LABEL PROPAGATION
##############################################################################

lp_tr = cluster_label_prop(jaccard_graph_tr, 
                         weights = jaccard_graph_tr$weights)

lp_knn = cluster_label_prop(jaccard_graph_knn, 
                         weights = jaccard_graph_knn$weights)
par(mfrow=c(1,2))
plotar(lp_tr, jaccard_graph_tr, "Label Propagation threshold")
plotar(lp_knn, jaccard_graph_knn, "Label Propagation knn")

lp_tr
membership(lp_tr)
sizes(lp_tr)
modularity(lp_tr)
communities(lp_tr)
is_hierarchical(lp_tr)


##############################################################################
# 4. WALKTRAP
##############################################################################

wt_tr = cluster_walktrap(jaccard_graph_tr, 
                         weights = jaccard_graph_tr$weights,
                         steps = 3,
                         merges = TRUE,
                         modularity = TRUE, 
                         membership = TRUE)

wt_knn = cluster_walktrap(jaccard_graph_knn, 
                          weights = jaccard_graph_knn$weights,
                          steps = 3,
                          merges = TRUE,
                          modularity = TRUE, 
                          membership = TRUE)
par(mfrow=c(1,2))
plotar(wt_tr, jaccard_graph_tr, "Walktrap threshold")
plotar(wt_knn, jaccard_graph_knn, "Walktrap knn")

wt_tr
membership(wt_tr)
sizes(wt_tr)
modularity(wt_tr)
communities(wt_tr)
is_hierarchical(wt_tr)
dend_wt = as.dendrogram(wt_tr)
hc_wt = as.hclust(wt_tr)
hc_wt$order
plot_dendrogram(wt_tr)
title(cex.main= 0.8, main = "Dendrogram Walktrap Threshold")


##############################################################################
# 5. LEADING EIGENVECTOR
##############################################################################

le_tr = cluster_leading_eigen(jaccard_graph_tr, 
                              weights = jaccard_graph_tr$weights)

le_knn = cluster_leading_eigen(jaccard_graph_knn, 
                               weights = jaccard_graph_knn$weights)
par(mfrow=c(1,2))
plotar(le_tr, jaccard_graph_tr, "Leading Eigenvector threshold")
plotar(le_knn, jaccard_graph_knn, "Leading Eigenvector knn")

le_tr
membership(le_tr)
sizes(le_tr)
modularity(le_tr)
communities(le_tr)
is_hierarchical(le_tr)
dend_le = as.dendrogram(le_tr)
hc_le = as.hclust(le_tr)
hc_le$order
plot_dendrogram(le_tr)
title(cex.main= 0.8, main = "Dendrogram Leading \n Eigenvector Treshold")



##############################################################################
# 6. Optimal
##############################################################################
op_tr = cluster_optimal(jaccard_graph_tr, 
                        weights = jaccard_graph_tr$weights)

op_knn = cluster_optimal(jaccard_graph_knn, 
                        weights = jaccard_graph_knn$weights)
par(mfrow=c(1,2))
plotar(op_tr, jaccard_graph_tr, "Optimal threshold")
plotar(op_knn, jaccard_graph_knn, "Optimal Eigenvector knn")

op_tr
membership(op_tr)
sizes(op_tr)
modularity(op_tr)
communities(op_tr)
is_hierarchical(op_tr)


##############################################################################
# 7. Louvain
##############################################################################
lo_tr = cluster_louvain(jaccard_graph_tr, 
                        weights = jaccard_graph_tr$weights)

lo_knn = cluster_louvain(jaccard_graph_knn, 
                        weights = jaccard_graph_knn$weights)
par(mfrow=c(1,2))
plotar(lo_tr, jaccard_graph_tr, "Louvain threshold")
plotar(lo_knn, jaccard_graph_knn, "Louvain Eigenvector knn")

lo_tr
membership(lo_tr)
sizes(lo_tr)
modularity(lo_tr)
communities(lo_tr)
is_hierarchical(lo_tr)


##############################################################################
# 8. Fast Greedy
##############################################################################
fg_tr = cluster_fast_greedy(jaccard_graph_tr, 
                            merges = TRUE,
                            modularity = TRUE,
                            membership = TRUE,
                            weights = jaccard_graph_tr$weights)

fg_knn = cluster_fast_greedy(jaccard_graph_knn, merges = TRUE,
                             modularity = TRUE,
                             membership = TRUE,
                             weights = jaccard_graph_knn$weights)

plotar(fg_tr, jaccard_graph_tr, "Fast Greedy threshold")
plotar(fg_knn, jaccard_graph_knn, "Fast Greedy knn")



##############################################################################
# GRAPH PROPERTIES                                                           #
##############################################################################

# The proportion of present edges from all possible edges in the network.
edge_density(jaccard_graph_tr, loops=F) 

# 
transitivity(jaccard_graph_tr, type="global")  

# The length of the shortest path between two nodes
diameter(jaccard_graph_tr, directed=F, weights=NA)

# 
degree(jaccard_graph_tr, mode="all")

# centrality based on distance to others in the graph
# Inverse of the node’s average geodesic distance to others in the network
closeness(jaccard_graph_tr, mode="all", weights=NA) 

# centrality proportional to the sum of connection centralities
# Values of the first eigenvector of the graph matrix
eigen_centrality(jaccard_graph_tr, directed=T, weights=NA)

# centrality based on a broker position connecting others
# Number of geodesics that pass through the node or the edge
betweenness(jaccard_graph_tr, directed=T, weights=NA)

# 
cliques(jaccard_graph_tr)

# Assortativity and Homophily: the tendency of nodes to connect to others 
# who are similar on some variable.
assortativity_degree(jaccard_graph_tr, directed=F)

