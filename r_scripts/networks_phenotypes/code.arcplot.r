#!/usr/bin/env Rscript


# creates arcdiagram

# first run bash script text_to_gml.sh
# ex:
# bash ./bash_scripts/text_to_gml.sh ./files/stomach_skyblue_males_mean_expression.txt ./files/stomach_skyblue_module_males_edges.txt ./files/stomach_skyblue_module_males_edges.gml


library(reshape)
library(plyr)
library(arcdiagram)
library(igraph)


#args = c("./files/stomach_skyblue_module_males_edges.gml", "./plots/wgcna_networks_traits/stomach_skyblue_module_males_edges.pdf")
#args = c("./files/stomach_skyblue_module_females_edges.gml", "./plots/wgcna_networks_traits/stomach_skyblue_module_females_edges.pdf")
args = c("./files/thyroid_greenyellow_module_males_edges.gml", "./plots/wgcna_networks_traits/thyroid_greenyellow_module_males_edges.pdf")
#args = c("./files/thyroid_greenyellow_module_females_edges.gml", "./plots/wgcna_networks_traits/thyroid_greenyellow_module_females_edges.pdf")
#args = commandArgs(trailingOnly=TRUE)


# location of 'gml' file
mis_file = args[1]

# read 'gml' file
mis_graph = read.graph(mis_file, format="gml")

# get edgelist
edgelist = get.edgelist(mis_graph)

# get vertex labels
vlabels = get.vertex.attribute(mis_graph, "label")

# get vertex groups
vgroups = get.vertex.attribute(mis_graph, "group")

# get vertex fill color
vfill = get.vertex.attribute(mis_graph, "fill")

# get vertex border color
vborders = get.vertex.attribute(mis_graph, "border")

# get vertex degree
degrees = degree(mis_graph)

# get edges value
values = get.edge.attribute(mis_graph, "value")

# data frame with vgroups, degree, vlabels and ind
x = data.frame(vgroups, degrees, vlabels, ind=1:vcount(mis_graph))

# arranging by vgroups and degrees
#y = arrange(x, vlabels)
y = arrange(x, degrees)
#y = arrange(x, desc(degrees))

# get ordering 'ind'
new_ord = y$ind

#png(args[2],w=2000,h=1500,res=70)
pdf(args[2],w=16,h=8)
par(oma=c(4,0,0,0))
arcplot(edgelist, vertices = get.vertex.attribute(mis_graph)$id, ordering=new_ord, labels=vlabels,  cex.labels=1.5,
  show.nodes=TRUE, col.nodes=vborders, bg.nodes="#de2d26", col.labels="black", cex.nodes = (log(degrees)+0.5)*2,
  pch.nodes=21, lwd.nodes = 1, line=-0.5, col.arcs = ifelse(values > 7, "#33333380", "#00000000"), lwd.arcs = 2)
dev.off()

#col.arcs = hsv(0, 0, 0.2, 0.5)
#col.arcs = gray(abs((values/10)-1))
#lwd.arcs = 0.2*values
#cex.nodes = 5
