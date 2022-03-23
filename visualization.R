library(igraph)
library(RColorBrewer)
help("brewer.pal")
edges <- read.csv("edge.csv", sep =",")
attrs <- read.csv("attr.csv", sep =",")

#sample_idx -> sample(nrow(edges), 10000)
#edges <- edges[sample_idx, ]
#attrs <- attrs[sample_idx, ]

actorGraph <- graph.data.frame(edges)
V(actorGraph)
vcount(actorGraph)
actorGraph <- set_vertex_attr(actorGraph, "group", value = attrs$group)
actorGraph <- set_vertex_attr(actorGraph, "salesrank", value = attrs$salesrank)
actorGraph <- set_vertex_attr(actorGraph, "reviews_num", value = attrs$reviews_num)
actorGraph <- set_vertex_attr(actorGraph, "avg_rating", value = attrs$avg_rating)

comp <- components(actorGraph)
comp
n<-length(comp$csize)
giantGraph <- actorGraph %>% 
  induced.subgraph(., which(comp$membership == which(comp$csize == sort(comp$csize)[n-2])))

col <- data.frame(group = unique(vertex_attr(giantGraph, "group")), stringsAsFactors = F)

if (length(col$group) >= 3) {
  col$color <- brewer.pal(nrow(col), "GnBu")
} else if (length(col$group) == 2){
  col$color <- c("Yellow","Green")
} else {
  col$color <- c("Yellow")
}

V(giantGraph)$color <- col$color[match(V(giantGraph)$group, col$group)]

giantGraph %>% 
  plot(.,
       layout = layout_with_kk(.), ## force-directed graph layout
       edge.arrow.size = .3,
       vertex.label = NA,
       vertex.size = 4,
       vertex.label.cex = .5,
       vertex.label.color = 'black')

giantGraph %>% 
  plot(.,
       layout = layout_with_kk(.),
       # layout = layout_with_sugiyama(.),
       edge.arrow.size = .3,
       vertex.size = 4,
       vertex.label = NA,
       vertex.color = adjustcolor(graph.coreness(.), alpha.f = .3),
       vertex.label.cex = .5,
       vertex.label.color = 'black',
       mark.groups = by(seq_along(graph.coreness(.)), graph.coreness(.), invisible),
       mark.shape = 1/4,
       mark.col = rainbow(length(unique(graph.coreness(.))),alpha = .1),
       mark.border = NA
  )

vcount(giantGraph)

actorGraph %>% degree.distribution(.,mode="in") %>% 
  plot(., col = 'black', pch = 19, cex = 1.5,
       main = 'In-degree Distribution',
       ylab = 'Density',
       xlab = 'In-degree')
actorGraph %>% 
  degree.distribution(.,cumulative = TRUE,mode ='in') %>% 
  plot(1:(max(degree(actorGraph,mode='in'))+1),., ## since log doesn't take 0, add 1 to every degree
       log='xy', type = 'l',
       main = 'Log-Log Plot of In-degree',
       ylab = 'CCDF',
       xlab = 'In-degree')

actorGraph %>% degree.distribution(.,mode="out") %>% 
  plot(., col = 'black', pch = 19, cex = 1.5,
       main = 'Out-degree Distribution',
       ylab = 'Density',
       xlab = 'Out-degree')
# Plot a log-log plot
actorGraph %>% 
  degree.distribution(.,cumulative = TRUE,mode ='out') %>% 
  plot(1:(max(degree(actorGraph,mode='out'))+1), ## since log doesn't take 0, add 1 to every degree
       ., log='xy', type = 'l',
       main = 'Log-Log Plot of Out-degree',
       ylab = 'CCDF',
       xlab = 'Out-degree')
