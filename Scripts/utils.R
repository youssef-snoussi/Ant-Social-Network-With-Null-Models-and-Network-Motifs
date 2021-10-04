########################################################################################
# In this file there are definition of some necessary functions that we will use in    #
# main program, for that we're going to need the following packages, if not installed  #
# incomment the following lines and run:                                               #
########################################################################################

####### Igraph packge #######
# install.packages("igraph")

###### Comprehenr package #######
# install.packages("comprehenr")



# Creating graphs from edges files
create_graphs <- function(path, directed=F, columns, sep=" "){
  
  library(comprehenr)
  dataframes <- to_list(for (f in list.files(path)) read.table(paste0(path,"/",f), col.names=columns, sep=sep))
  
  library(igraph)
  graphList <- to_list(for (df in dataframes) graph_from_data_frame(df, directed=directed))
  
  return(graphList)
  }

# Graph stats
graph_stats <- function(graphList){
  gnodes <- to_vec(for (g in graphList) vcount(g))               # Total nodes per graph
  gedges <- to_vec(for (g in graphList) ecount(g))               # Totale edges prt graph
  gAD <- gedges / gnodes                                           # Average degree per graph
  gEW <- to_vec(for (g in graphList) mean(E(g)$weight))          # Edge weights
  gdeg <- to_list(for (g in graphList) degree(g))                # Node degrees
  gAPL <- to_vec(for (g in graphList) average.path.length(g))    # Average path length per graph
  gDiam <- to_vec(for (g in graphList) diameter(g))              # Graphs' diameters
  gDs <- to_vec(for (g in graphList) edge_density(g))            # Edge Density
  
  graphs_summary <- as.data.frame(cbind(gAD,gEW,gAPL,gDs))
  colnames(graphs_summary) <- c("AverageDegree", "EdgeWeight", "AveragePathLength", "Density")
  return(graphs_summary)
}

# Creating dataframes containing graphs' informations like degree etc
create_dataframe_from_graph <- function(graph){
  library(igraph)
  edges <- to_list(for (e in E(graph)) ends(graph,e))
  ed_wgt <- edge.attributes(graph)[["weight"]]
  ed_str <- to_vec(for (e in edges) (strength(graph,v=as.character(e[1])) + strength(graph,as.vector(e[2])))/2)
  ed_bet <- edge_betweenness(graph)
  edges <- get.edgelist(graph)
  df <- as.data.frame(cbind(edges, ed_wgt, ed_str, ed_bet), col.names = c("V1","V2","weight","degree","betweeness"))
  df$ed_str <- as.numeric(as.character(df$ed_str))
  df$ed_bet <- as.numeric(as.character(df$ed_bet))
  df$ed_wgt <- as.numeric(as.character(df$ed_wgt))
  return(df)
}

# Creating nul  l models with the model of Erdos-Renyi random model
nulmodel <- function(graph, weighted=T, directed=F){
  library(igraph)
  
  # Number of vertices
  nv <- vcount(graph)
  
  # Number of edges
  ne <- ecount(graph) 
  
  # Grab the weights and permute them
  wgt <- edge.attributes(graph)[["weight"]]
  wgt <- sample(max(wgt), size=length(wgt), replace=T)
  
  # Create a random graph with Erdos-Renyi model
  g <- sample_gnm(nv, ne, directed=directed)
  g <- set.edge.attribute(g, name="weight", value=wgt)
  V(g)$name <- as.character(1:vcount(g))
  return(g)
}

# search triad motifs
triad_cli <- function(graph){
  library(igraph)
  cl.tri <- cliques(graph, min=3, max=3)
  df <- lapply(cl.tri, function(x){V(graph)$name[x]})
  return(df)
}

#calculating subgraphs' weights defined as the sum of weights of each edge
sub_weights <- function(graph, tr){
  wgt <- edge.attributes(graph)[["weight"]]
  ed1 <- get.edge.ids(graph, c(tr[1],tr[2]))
  ed2 <- get.edge.ids(graph, c(tr[2],tr[3]))
  ed3 <- get.edge.ids(graph, c(tr[3],tr[1]))
  sub_wgt <- wgt[ed1]+wgt[ed2]+wgt[ed3]
}





