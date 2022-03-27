##                                                      ##
##    Source:
##    Katya Ognyanova, katya@ognyanova.net              ##
##    www.kateto.net/netscix2016                        ##
##                                                      ##
##======================================================##



# slide_p2

g1 <- graph(edges=c(1,2, 2,1,2,3, 3,2,3,1,4,6,6,4,3,7,2,5,2,6), n=7, directed=T )
plot(g1,vertex.color='skyblue',vertex.label.family='roboto',
     edge.arrow.size=0.3)

# ================ 2. Networks in igraph ================

rm(list = ls()) # Remove all the objects we created so far.

library(igraph) # Load the igraph package

g <- graph( c("John", "Jim", "Jim", "Jack", "Jim", "Jack", "John", "John"), 
             isolates=c("Jesse", "Janis", "Jennifer", "Justin") )  
# In named graphs we can specify isolates by providing a list of their names.

plot(g, edge.arrow.size=.5, vertex.color="gold", vertex.size=15, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=1, edge.curved=0.2) 


#  ------->> Edge, vertex, and network attributes --------

# Access vertices and edges:
E(g) # The edges of the object
V(g) # The vertices of the object

# You can also manipulate the network matrix directly:
g[]
g[1,]
g[3,3] <- 10
g[5,7] <- 10
g[]

# Add attributes to the network, vertices, or edges:
V(g)$name # automatically generated when we created the network.
V(g)$gender
V(g)$gender <- c("male", "male", "male", "male", "female", "female", "male")
E(g)$type <- "email" # Edge attribute, assign "email" to all edges
E(g)$weight <- 10    # Edge weight, setting all existing edges to 10

# Let's take a look at the description of the igraph object.
# Those will typically start with up to four letters:
# 1. D or U, for a directed or undirected graph
# 2. N for a named graph (where nodes have a name attribute)
# 3. W for a weighted graph (where edges have a weight attribute)
# 4. B for a bipartite (two-mode) graph (where nodes have a type attribute)
#
# The two numbers that follow refer to the number of nodes and edges in the graph. 
# The description also lists graph, node & edge attributes, for example:
# (g/c) - graph-level character attribute
# (v/c) - vertex-level character attribute
# (e/n) - edge-level numeric attribute

# Delete attributes
g
g <- delete_edge_attr(g,"weight")
g <- delete_vertex_attr(g,"gender")
g

plot(g, edge.arrow.size=.1, vertex.label.color="black", vertex.label.family='font',
     vertex.color=c("pink", "skyblue")[1+(V(g)$gender=="male")] ) 

# g has two edges going from Jim to Jack, and a loop from John to himself.
# We can simplify our graph to remove loops & multiple edges between the same nodes.

?simplify
gs <- simplify( g)

# EXERCISE
# why doesn't this code change vertex size??
plot(gs, edge.arrow.size=.1, vertex.label.color="black", vertex.lable.cex=0.8,
     vertex.label.family='font',
     vertex.color=c("pink", "skyblue")[1+(V(g)$gender=="male")])
gs



# ------->> Specific graphs and graph models --------

# Empty graph
eg <- make_empty_graph(40)
plot(eg, vertex.size=10, vertex.label=NA)

# Full graph
fg <- make_full_graph(40)
plot(fg, vertex.size=10, vertex.label=NA)

# Star graph 
st <- make_star(40)
plot(st, vertex.size=10, vertex.label=NA) 

# Tree graph
tr <- make_tree(40, children = 3, mode = "undirected")
plot(tr, vertex.size=10, vertex.label=NA) 

# Ring graph
rn <- make_ring(40)
plot(rn, vertex.size=10, vertex.label=NA)

# Erdos-Renyi random graph 
# ('n' is number of nodes, 'm' is the number of edges)
er <- sample_gnm(n=100, m=40) 
plot(er, vertex.size=6, vertex.label=NA)  

# Watts-Strogatz small-world graph
# Creates a lattice with 'dim' dimensions of 'size' nodes each, and rewires edges 
# randomly with probability 'p'. You can allow 'loops' and 'multiple' edges.
# The neighborhood in which edges are connected is 'nei'.
sw <- sample_smallworld(dim=2, size=10, nei=1, p=0.1)
plot(sw, vertex.size=6, vertex.label=NA, layout=layout_in_circle)
 
# Barabasi-Albert preferential attachment model for scale-free graphs
# 'n' is number of nodes, 'power' is the power of attachment (1 is linear)
# 'm' is the number of edges added on each time step 
 ba <-  sample_pa(n=100, power=1, m=1,  directed=F)
 plot(ba, vertex.size=6, vertex.label=NA)
 
#igraph can also give you some notable historical graphs. For instance:
 zach <- graph("Zachary") # the Zachary carate club
 plot(zach, vertex.size=10, vertex.label=NA)
 
  # Rewiring a graph
 # 'each_edge()' is a rewiring method that changes the edge endpoints
 # uniformly randomly with a probability 'prob'.
 rn.rewired <- rewire(rn, each_edge(prob=0.1))
 plot(rn.rewired, vertex.size=10, vertex.label=NA)
 
 # Rewire to connect vertices to other vertices at a certain distance. 
 rn.neigh = connect.neighborhood(rn, 5)
 plot(rn.neigh, vertex.size=8, vertex.label=NA) 
 
 
 # Combine graphs (disjoint union, assuming separate vertex sets): %du%
 plot(rn, vertex.size=10, vertex.label=NA) 
 plot(tr, vertex.size=10, vertex.label=NA) 
 plot(rn %du% tr, vertex.size=10, vertex.label=NA) 

  
 
# ================ 3. Reading network data from files ================

 
rm(list = ls()) # clear the workspace again

# Download the archive with the data files from http://bitly.com/netscix2016 
 
# Set the working directory to the folder containing the workshop files:
setwd("C:/DOCS/Conferences/2016-NetSciX/NetSciX Workshop")  
 
# DATASET 1: edgelist 

nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)

# Examine the data:
head(nodes)
head(links)
nrow(nodes); length(unique(nodes$id))
nrow(links); nrow(unique(links[,c("from", "to")]))

# Collapse multiple links of the same type between the same two nodes
# by summing their weights, using aggregate() by "from", "to", & "type":
# (we don't use "simplify()" here so as not to collapse different link types)
links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL


# DATASET 2: matrix 

nodes2 <- read.csv("Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
links2 <- read.csv("Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)

# Examine the data:
head(nodes2)
head(links2)

# links2 is an adjacency matrix for a two-mode network:
links2 <- as.matrix(links2)
dim(links2)
dim(nodes2)


# ================ 4. Turning networks into igraph objects ================ 
 
library(igraph)

#  ------->> DATASET 1 -------- 

# Converting the data to an igraph object:
# The graph.data.frame function, which takes two data frames: 'd' and 'vertices'.
# 'd' describes the edges of the network - it should start with two columns 
# containing the source and target node IDs for each network tie.
# 'vertices' should start with a column of node IDs.
# Any additional columns in either data frame are interpreted as attributes.

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 

# Examine the resulting object:
class(net)
net 

# We can look at the nodes, edges, and their attributes:
E(net)
V(net)
E(net)$type
V(net)$media

plot(net, edge.arrow.size=.2,vertex.label=NA)

# Removing loops from the graph:
net <- simplify(net, remove.multiple = T, remove.loops = T) 

# If you need them, you can extract an edge list or a matrix from igraph networks.
as_edgelist(net, names=T)
as_adjacency_matrix(net, attr="weight")

# Or data frames describing nodes and edges:
as_data_frame(net, what="edges")
as_data_frame(net, what="vertices")

# ================ 6. Network and node descriptives ================


# Density
# The proportion of present edges from all possible ties.
edge_density(net, loops=F)
ecount(net)/(vcount(net)*(vcount(net)-1)) #for a directed network

# Reciprocity
# The proportion of reciprocated ties (for a directed network).
reciprocity(net)
dyad_census(net) # Mutual, asymmetric, and null node pairs
2*dyad_census(net)$mut/ecount(net) # Calculating reciprocity

# Transitivity
# global - ratio of triangles (direction disregarded) to connected triples
# local - ratio of triangles to connected triples each vertex is part of
transitivity(net, type="global")  # net is treated as an undirected network
transitivity(as.undirected(net, mode="collapse")) # same as above
transitivity(net, type="local")
triad_census(net) # for directed networks

# Triad types (per Davis & Leinhardt):
# 
# 003  A, B, C, empty triad.
# 012  A->B, C 
# 102  A<->B, C  
# 021D A<-B->C 
# 021U A->B<-C 
# 021C A->B->C
# 111D A<->B<-C
# 111U A<->B->C
# 030T A->B<-C, A->C
# 030C A<-B<-C, A->C.
# 201  A<->B<->C.
# 120D A<-B->C, A<->C.
# 120U A->B<-C, A<->C.
# 120C A->B->C, A<->C.
# 210  A->B<->C, A<->C.
# 300  A<->B<->C, A<->C, completely connected.


# Diameter (longest geodesic distance)
# Note that edge weights are used by default, unless set to NA.
diameter(net, directed=F, weights=NA)
diameter(net, directed=F)
diam <- get_diameter(net, directed=T)
diam

# Note: vertex sequences asked to behave as a vector produce numeric index of nodes
class(diam)
as.vector(diam)

# Color nodes along the diameter:
vcol <- rep("gray40", vcount(net))
vcol[diam] <- "gold"
ecol <- rep("gray80", ecount(net))
ecol[E(net, path=diam)] <- "orange" 
# E(net, path=diam) finds edges along a path, here 'diam'
plot(net, vertex.color=vcol, edge.color=ecol, edge.arrow.mode=0)

# Node degrees
# 'degree' has a mode of 'in' for in-degree, 'out' for out-degree,
# and 'all' or 'total' for total degree. 
deg <- degree(net, mode="all")
plot(net, vertex.size=deg*3)
hist(deg, breaks=1:vcount(net)-1, main="Histogram of node degree")

# Degree distribution
deg.dist <- degree_distribution(net, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")


# Centrality & centralization

# Centrality functions (vertex level) and centralization functions (graph level).
# The centralization functions return "res" - vertex centrality, "centralization", 
# and "theoretical_max" - maximum centralization score for a graph of that size.
# The centrality functions can run on a subset of nodes (set with the "vids" parameter)

# Degree (number of ties)
degree(net, mode="in")
centr_degree(net, mode="in", normalized=T)

# Closeness (centrality based on distance to others in the graph)
# Inverse of the node's average geodesic distance to others in the network
closeness(net, mode="all", weights=NA) 
centr_clo(net, mode="all", normalized=T) 

# Eigenvector (centrality proportional to the sum of connection centralities)
# Values of the first eigenvector of the graph adjacency matrix
eigen_centrality(net, directed=T, weights=NA)
centr_eigen(net, directed=T, normalized=T) 

# Betweenness (centrality based on a broker position connecting others)
# (Number of geodesics that pass through the node or the edge)
betweenness(net, directed=T, weights=NA)
edge_betweenness(net, directed=T, weights=NA)
centr_betw(net, directed=T, normalized=T)



-----------------------------------
# * TASK FOR WORKSHOP PARTICIPANTS:

# Compute the degree, closeness, eigenvector, and betweenness centrality of
# the actors in the Zachary karate club network. Plot the network, sizing the
# nodes based on the different centrality types.

-----------------------------------

  

# Hubs and authorities

# The hubs and authorities algorithm developed by Jon Kleinberg was initially used 
# to examine web pages. Hubs were expected to contain catalogues with a large number 
# of outgoing links; while authorities would get many incoming links from hubs, 
# presumably because of their high-quality relevant information. 

hs <- hub_score(net, weights=NA)$vector
as <- authority_score(net, weights=NA)$vector

par(mfrow=c(1,2))
 plot(net, vertex.size=hs*50, main="Hubs")
 plot(net, vertex.size=as*30, main="Authorities")
dev.off()


# ================ 7. Distances and paths ================


# Average path length 
# The mean of the shortest distance between each pair of nodes in the network 
# (in both directions for directed graphs). 
mean_distance(net, directed=F)
mean_distance(net, directed=T)

# We can also find the length of all shortest paths in the graph:
distances(net) # with edge weights
distances(net, weights=NA) # ignore weights

# We can extract the distances to a node or set of nodes we are interested in.
# Here we will get the distance of every media from the New York Times.
dist.from.NYT <- distances(net, v=V(net)[media=="NY Times"], to=V(net), weights=NA)

# Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.NYT)+1)
col <- col[dist.from.NYT+1]

plot(net, vertex.color=col, vertex.label=dist.from.NYT, edge.arrow.size=.6, 
     vertex.label.color="white")

# We can also find the shortest path between specific nodes.
# Say here between MSNBC and the New York Post:
news.path <- shortest_paths(net, 
                            from = V(net)[media=="MSNBC"], 
                             to  = V(net)[media=="New York Post"],
                             output = "both") # both path nodes and edges

# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(net))
ecol[unlist(news.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(net))
ew[unlist(news.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(net))
vcol[unlist(news.path$vpath)] <- "gold"

plot(net, vertex.color=vcol, edge.color=ecol, 
     edge.width=ew, edge.arrow.mode=0)


# Identify the edges going into or out of a vertex, for instance the WSJ.
# For a single node, use 'incident()', for multiple nodes use 'incident_edges()'
inc.edges <- incident(net, V(net)[media=="Wall Street Journal"], mode="all")

# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(net))
ecol[inc.edges] <- "orange"
vcol <- rep("grey40", vcount(net))
vcol[V(net)$media=="Wall Street Journal"] <- "gold"
plot(net, vertex.color=vcol, edge.color=ecol)


# We can also easily identify the immediate neighbors of a vertex, say WSJ.
# The 'neighbors' function finds all nodes one step out from the focal actor.
# To find the neighbors for multiple nodes, use 'adjacent_vertices()'.
# To find node neighborhoods going more than one step out, use function 'ego()'
# with parameter 'order' set to the number of steps out to go from the focal node(s).

neigh.nodes <- neighbors(net, V(net)[media=="Wall Street Journal"], mode="out")

# Set colors to plot the neighbors:
vcol[neigh.nodes] <- "#ff9d00"
plot(net, vertex.color=vcol)

# Special operators for the indexing of edge sequences: %--%, %->%, %<-%
# E(network)[X %--% Y] selects edges between vertex sets X and Y, ignoring direction
# E(network)[X %->% Y] selects edges from vertex sets X to vertex set Y
# E(network)[X %->% Y] selects edges from vertex sets Y to vertex set X

# For example, select edges from newspapers to online sources:
E(net)[ V(net)[type.label=="Newspaper"] %->% V(net)[type.label=="Online"] ]

# Cocitation (for a couple of nodes, how many shared nominations they have)
cocitation(net)



# ================ 8. Subgroups and communities ================

# Converting 'net' to an undirected network.
# There are several ways to do that: we can create an undirected link between any pair
# of connected nodes (mode="collapse), or create an undirected link for each directed
# one (mode="each"), or create an undirected link for each symmetric link (mode="mutual").
# In cases when A -> B and B -> A are collapsed into a single undirected link, we
# need to specify what to do with the edge attributes. Here we have said that
# the 'weight' of links should be summed, and all other edge attributes ignored.

net.sym <- as.undirected(net, mode="collapse", edge.attr.comb=list(weight="sum", "ignore"))


#  ------->> Cliques --------

 # Find cliques (complete subgraphs of an undirected graph)
cliques(net.sym) # list of cliques       
sapply(cliques(net.sym), length) # clique sizes
largest_cliques(net.sym) # cliques with max number of nodes

vcol <- rep("grey80", vcount(net.sym))
vcol[unlist(largest_cliques(net.sym))] <- "gold"
plot(net.sym, vertex.label=V(net.sym)$name, vertex.color=vcol)



#  ------->> Communities --------

# A number of algorithms aim to detect groups that consist of densely connected nodes
# with fewer connections across groups. 

# Community detection based on edge betweenness (Newman-Girvan)
# High-betweenness edges are removed sequentially (recalculating at each step)
# and the best partitioning of the network is selected.
ceb <- cluster_edge_betweenness(net) 
dendPlot(ceb, mode="hclust")
plot(ceb, net) 

# Let's examine the community detection igraph object:
class(ceb)
length(ceb)     # number of communities
membership(ceb) # community membership for each node
crossing(ceb, net)   # boolean vector: TRUE for edges across communities
modularity(ceb) # how modular the graph partitioning is

# High modularity for a partitioning reflects dense connections within communities 
# and sparse connections across communities.


# Community detection based on propagating labels
# Assigns node labels, randomizes, and replaces each vertex's label with
# the label that appears most frequently among neighbors. Repeated until
# each vertex has the most common label of its neighbors.
clp <- cluster_label_prop(net)
plot(clp, net)

# Community detection based on greedy optimization of modularity
cfg <- cluster_fast_greedy(as.undirected(net))
plot(cfg, as.undirected(net))
 
# We can also plot the communities without relying on their built-in plot:
V(net)$community <- cfg$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(net, vertex.color=colrs[V(net)$community])


-----------------------------------
# * TASK FOR WORKSHOP PARTICIPANTS:

# Plot the results of the three different community detection algorithms 
#  applied to  the Zachary karate club network. 
  
-----------------------------------


# K-core decomposition
# The k-core is the maximal subgraph in which every node has degree of at least k
# This also means that the (k+1)-core will be a subgraph of the k-core.
# The result here gives the coreness of each vertex in the network.
kc <- coreness(net, mode="all")
plot(net, vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc])


# ================ 9. Assortativity and Homophily ================

# Assortativity (homophily)
# The tendency of nodes to connect to others who are similar on some variable.
# assortativity_nominal() is for categorical variables (labels)
# assortativity() is for ordinal and above variables
# assortativity_degree() checks assortativity in node degrees

V(net)$type.label
V(net)$media.type

assortativity_nominal(net, V(net)$media.type, directed=F)

assortativity(net, V(net)$audience.size, directed=F)

assortativity_degree(net, directed=F)


# ================ The End ================




