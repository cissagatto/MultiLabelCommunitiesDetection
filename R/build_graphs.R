##################################################################################################
# Tests Communities Detection Methods                                                            #
# Copyright (C) 2021                                                                             #
#                                                                                                #
# This code is free software: you can redistribute it and/or modify it under the terms of the    #
# GNU General Public License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version. This code is distributed in the hope       #
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    #
# more details.                                                                                  #
#                                                                                                #
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin                     #
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus Sao Carlos           #
# Computer Department (DC: https://site.dc.ufscar.br/)                                           #
# Program of Post Graduation in Computer Science (PPG-CC: http://ppgcc.dc.ufscar.br/)            #
# Bioinformatics and Machine Learning Group (BIOMAL: http://www.biomal.ufscar.br/)               #
#                                                                                                #
##################################################################################################

rm(list = ls()) # Remove all the objects we created so far.

##################################################################################################
# Configures the workspace according to the operating system                                     #
##################################################################################################
sistema = c(Sys.info())
FolderRoot = ""
if (sistema[1] == "Linux"){
  FolderRoot = paste("/home/", sistema[7], "/TestsCommunitiesDetection", sep="")
  setwd(FolderRoot)
} else {
  FolderRoot = paste("C:/Users/", sistema[7], "/TestsCommunitiesDetection", sep="")
  setwd(FolderRoot)
}
setwd(FolderRoot)

##################################################################################################
FolderScripts =  paste(FolderRoot, "/R", sep="")

setwd(FolderScripts)
source("libraries.R")

setwd(FolderScripts)
source("utils.R")


##################################################################################################
# Options Configuration                                                                          #
##################################################################################################
options(java.parameters = "-Xmx32g")
options(show.error.messages = TRUE)
options(scipen=10)


##################################################################################################
# Read the dataset file with the information for each dataset                                    #
##################################################################################################
setwd(FolderRoot)
datasets <- data.frame(read.csv("datasets.csv"))
ds <- datasets[29,]
dataset_name = ds$Name


##################################################################################################
# Open dataset
##################################################################################################
folder = folders(dataset_name)
setwd(folder$FolderTR)
dados_1 = foreign::read.arff("GpositiveGO-Split-Tr-1.arff")

# separando rótulos e instancias
labels = dados_1[,ds$LabelStart:ds$LabelEnd]
instances = dados_1[,ds$AttStart:ds$AttEnd]

#obtendo nomes dos rótulos e instâncias
setwd(folder$FolderNamesLabels)
nomes_labels = data.frame(read.csv(paste(dataset_name, "-NamesLabels.csv", sep="")))
colnames(nomes_labels) = c("i", "nomes")

# abrindo os resultados das medidas de similaridades
setwd(folder$FolderSimilarities)
sim = read.csv("vari.csv", stringsAsFactors = FALSE)
rownames(sim) = sim$X
sim = sim[,-1]
sim_2 = as.matrix(sim)

# reorganizando as arestas e vértices
arestas = data.frame(pivot_longer(sim, cols = 1:4))
colnames(arestas) = c("to", "similarity")

from = c("Label1", "Label1", "Label1", "Label1",
         "Label2", "Label2", "Label2", "Label2",
         "Label3", "Label3", "Label3", "Label3",
         "Label4", "Label4", "Label4", "Label4")
arestas_2 = cbind(from, arestas)

#############################################################
# pesando as arestas com as probabilidades condicionais
setwd(folder$FolderSimilarities)
pesos = data.frame(read.csv("GpositiveGO-conditional-probabilities.csv"))
rownames(pesos) = pesos$X
pesos_2 = pesos[,-1]
pesos_3 = data.frame(pivot_longer(pesos_2, cols = 1:4))
weights = pesos_3$value
sim_3 = cbind(arestas_2, weights)

#############################################################
# pesando as arestas com os totais de cada classe
setwd(folder$FolderSimilarities)
pesos_2 = data.frame(read.csv("GpositiveGO-total-labels.csv"))
colnames(pesos_2) = c("", nomes_labels$nomes)
rownames(pesos_2) = nomes_labels$nomes
pesos_2 = pesos_2[,-1]
pesos_2 = data.frame(pivot_longer(pesos_2, cols = 1:4))
weights = pesos_2$value
sim_4 = cbind(arestas_2, weights)

###########################################
# EXEMPLO 1: grafos a partir de sim_2
grafo_sem_peso = graph.adjacency(sim_2, mode="undirected")
plot(grafo_sem_peso)

grafo_com_peso = graph.adjacency(sim_2, mode="undirected", weighted = TRUE)
plot(grafo_com_peso)

###########################################
# as_data_frame {igraph}	
# Creating igraph graphs from data frames or vice-versa
# This function creates an igraph graph from one or two data 
# frames containing the (symbolic) edge list and edge/vertex attributes.
# If vertices is NULL, then the first two columns of d are used as 
# a symbolic edge list and additional columns as edge attributes. 
# The names of the attributes are taken from the names of the columns.

# EXEMPLO 2: grafos a partir de sim_3
grafo_sim_3 = graph_from_data_frame(sim_3, directed=F) 
plot(grafo_sim_3)
E(grafo_sim_3) # arestas
V(grafo_sim_3) # vertices


###########################################
# EXEMPLO 2: grafos a partir de sim_4
grafo_sim_4 = graph_from_data_frame(sim_4, directed=FALSE) 
plot(grafo_sim_4)



# package ccccd
# Nearest Neighbor Graphs
# nearest neighbor, k-nearest neighbor, and mutual k-nearest neighbor (di)graphs.
# nng(x, dx, k, mutual, method, use.fnn, algorithm)
nng_1 = cccd::nng(x = sim_2, k=10, mutual = FALSE)


# rng {cccd}	R Documentation
# Relative Neighborhood Graph.
# the relative neighborhood graph defined by a set of points.
rng_1 = cccd::rng(x = sim_2, r = 2, k = 3)

# Average nearest neighbor degree
igraph::knn()

############################
# ego_size {igraph}	R Documentation
# Neighborhood of graph vertices
# These functions find the vertices not farther than a given limit 
# from another fixed vertex, these are called the neighborhood of the vertex.

# calculates the neighborhoods of the given vertices with the given order parameter.
grafo_ego = igraph::ego(grafo_com_peso, order=2, nodes = V(grafo_com_peso))

igraph::ego(sim_3, order=2, nodes = V(grafo_com_peso))

# calculates the size of the neighborhoods for the given vertices with the given order.
grafo_ego_size = igraph::ego_size(grafo_com_peso, order=2, nodes = V(grafo_com_peso))
plot(grafo_ego_size)

# s creates (sub)graphs from all neighborhoods of the given vertices 
# with the given order parameter. This function preserves the vertex,
# edge and graph attributes.
grafo_make_ego = igraph::make_ego_graph(grafo_com_peso, order= 2, nodes = V(grafo_com_peso))

grafo_ns = igraph::neighborhood.size(grafo_com_peso, order = 2,nodes = V(grafo_com_peso))
plot(grafo_ns)

grafo_n = igraph::neighborhood(grafo_com_peso, order = 2,nodes = V(grafo_com_peso))

##########################################
# neighbors {igraph}	R Documentation
# Neighboring (adjacent) vertices in a graph
# A vertex is a neighbor of another one (in other words, the two vertices 
# are adjacent), if they are incident to the same edge.
grafo_nn = igraph::neighbors(grafo_com_peso, v = V(grafo_com_peso))
plot(grafo_nn)

#############################################################
# COMUNIDADES
############################################################

# Spinglass
# weights = The weights of the edges. Either a numeric vector or NULL.
# If it is null and the input graph has a ‘weight’ edge attribute then 
# that will be used. If NULL and no such attribute is present then the 
# edges will have equal weights. Set this to NA if the graph was a ‘weight’ 
# edge attribute, but you don't want to use it for community detection. 
# A larger edge weight means a stronger connection for this function.
sc0 = spinglass.community(grafo_sim_3, implementation = "orig", update.rule = "simple")




dir = dir(folder$FolderSimilarities)
n = length(dir)
i = 1
for(i in 1:n){
  cat(i, " \t measure \t", dir[i])
  setwd(folder$FolderSimilarities)
  sim = read.csv(dir[i])
  rownames(sim) = sim$X
  sim = sim[,-1]
  sim_2 = as.matrix(sim)
  
  # grafo
  grafo1 = graph.adjacency(sim_2, mode="undirected")
  
  setwd(folder$FolderPlots)
  png(paste(dir[i], ".png", sep=""), width=1000, height=1000, res=200)
  print(plot(grafo1))
  dev.off()
  cat("\n")
  
  i = i + 1
  gc()
}