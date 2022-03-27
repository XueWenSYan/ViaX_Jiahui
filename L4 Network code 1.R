# This script heavily borrows from Dr. Anand Sokhey's ICPSR 2021 course examples.   

rm(list =ls())
install.packages("igraph")
install.packages("statnet")
install.packages("network")
library(igraph)
library(statnet,quietly = T)
library(network,quietly = T)
library(tidyverse)

m_game <- read.csv("C:\\Users\\xwuey\\Desktop\\relation.csv")
m_game[is.na(m_game)] <- 0
m_game<- m_game[,-1] %>% as.matrix() 
net_game <- network(m_game, matrix.type="adjacency");net_game
network.vertex.names(net_game) <- c("??????", "??????", "?????????", "D", "E")
gplot(net_game,displaylabels=T)

#create a network object from an adjacency matrix
netmat1 <- rbind(c(0,1,1,0,0),
                 c(0,0,1,1,0),
                 c(0,1,0,0,0),
                 c(0,0,0,0,0),
                 c(0,0,1,0,0))
rownames(netmat1) <- c("A", "B", "C", "D", "E")
(colnames(netmat1) <- c("A", "B", "C", "D", "E"));netmat1 
as.network(netmat1)
net1 <- network(netmat1, matrix.type="adjacency");net1

gplot(net1, vertex.col = 2, displaylabels = TRUE,vertex.cex = 0.2)


#create a network object from an edge list
netmat2 <- rbind(c(1,2),
                 c(1,3),
                 c(2,3),
                 c(2,4),
                 c(3,2),
                 c(5,3))
net2 <- network(netmat2, matrix.type = "edgelist")
network.vertex.names(net2) <- c("A", "B", "C", "D", "E")
summary(net2)
gplot(net2)
net2

#create graph objects from adjacency matrix or edge list
g1 <- igraph::graph_from_adjacency_matrix(netmat1) #calling netmat1 from above
g2 <- igraph::graph.edgelist(netmat2)
g1
plot(g2)

#plotting configurations
set.seed(1993)
plot(g1,vertex.size=30,vertex.color = 'coral3',
     vertex.frame.color=NA,vertex.label.family='arial',
     edge.arrow.size=0.4,layout=igraph::layout_with_fr(g1)) # https://igraph.org/r/doc/layout_.html; https://igraph.org/r/doc/plot.common.html

## more configs
?degree
plot(g1,vertex.size=sqrt(igraph::degree(g1,mode='out'))*25,vertex.label=NA,
     vertex.color = 'SkyBlue2',
     vertex.frame.color=NA,
     edge.arrow.size=0.2,layout=igraph::layout_as_star(g1)) # https://igraph.org/r

#importing and exporting an edgelist as data.frame
netmat3 <- rbind(c("A", "B"),
                 c("A", "C"),
                 c("B", "C"),
                 c("B", "D"),
                 c("C", "B"),
                 c("E", "C"))
class(netmat3)
typeof(netmat3)
net.df <- data.frame(netmat3)
net.df

#saving as external csv file
getwd()
write.csv(net.df, file = "C:/Users/xwuey/Desktop/MyData.csv")
setwd("G:/My Drive/ViaX/徐佳慧")

#importing csv file to R
net.edge <- read.csv(file = "C:/Users/xwuey/Desktop/MyData.csv")
net_import <- network(net.edge[,2:3], matrix.type = "edgelist")
summary(net_import)
