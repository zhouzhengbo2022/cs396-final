library(xtable)
library(mvtnorm)
library(igraph)
source('MultivarALAAMalt.R')
help("brewer.pal")
edges <- read.csv("edge.csv", sep =",")
attrs <- read.csv("attr.csv", sep =",")

#sample_idx -> sample(nrow(edges), 10000)
#edges <- edges[10000, ]
#attrs <- attrs[10000, ]

actorGraph <- graph.data.frame(edges)
V(actorGraph)
vcount(actorGraph)
actorGraph <- set_vertex_attr(actorGraph, "group", value = attrs$group)
actorGraph <- set_vertex_attr(actorGraph, "salesrank", value = attrs$salesrank)
actorGraph <- set_vertex_attr(actorGraph, "reviews_num", value = attrs$reviews_num)
actorGraph <- set_vertex_attr(actorGraph, "avg_rating", value = attrs$avg_rating)
actorGraph <- set_vertex_attr(actorGraph, "high_rating", value = attrs$high_rating)
actorGraph <- set_vertex_attr(actorGraph, "ranking10", value = attrs$ranking10)
actorGraph <- set_vertex_attr(actorGraph, "ranking20", value = attrs$ranking20)
actorGraph <- set_vertex_attr(actorGraph, "ranking50", value = attrs$ranking50)
actorGraph
comp <- components(actorGraph)
comp
n<-length(comp$csize)
giantGraph <- actorGraph %>% 
  induced.subgraph(., which(comp$membership == which(comp$csize == sort(comp$csize)[n-3])))

vcount(giantGraph)

data <- as_data_frame(giantGraph, "both")

adj <- as_adj(giantGraph)
adj <- as.matrix(adj)
n <- nrow(adj)

out.degree <-matrix( rowSums(adj), n, 1) # number of ties sent
in.degree <- matrix( colSums(adj) , n, 1 ) # number of ties received
rec.ties <-  matrix( rowSums(adj * t(adj) ), n , 1) # number of ties that are mutual
in.two.star <- matrix( choose(in.degree,2),n,1) #  in-stars refecting dispersion in popularity
out.two.star <- matrix( choose(out.degree,2),n,1) #  out-stars refecting dispersion in activity
mix.two.star <- in.degree*out.degree - rec.ties # correlation between indegree and outdegree
in.three.star <- matrix( choose(in.degree,3),n,1) # furhter measure of in-degree heterogeneity
out.three.star <- matrix( choose(out.degree,3),n,1) # furhter measure of out-degree heterogeneity
triangles <- rowSums( adj* (adj %*% t(adj) )  ) # embedded in transitive triads

covs <- cbind(out.degree, 
              in.degree,
              rec.ties,
              in.two.star,
              out.two.star,
              mix.two.star,
              in.three.star,
              out.three.star,
              triangles,
              data$vertices$reviews_num,
              data$vertices$avg_rating)
colnames(covs) <- c("indegree",
                    "outdegree",
                    "reciprochation" ,
                    "instar",
                    "outstar",
                    "twopath",
                    "in3star",
                    "out3star",
                    "transitive",
                    "reviews_num",
                    "avg_rating")

res.0 <- BayesALAAM(y = data$vertices$ranking10,           # dependent variable
                    ADJ = adj,           # network
                    covariates = covs,   # covariates
                    directed = TRUE,     # directed / undirecred network
                    Iterations = 10000,   # number of iterations
                    saveFreq = 100,      # print and save frequency
                    #contagion = 'none'
                    )  # type of contagion
plot(ts(res.0$Thetas))
plot(ts(res.0$Theta[,1:6]))
plot(ts(res.0$Theta[,7:11]))
## Summarize the results
write.res.table(burnin=1, # should be set sufficiently high
                datamat=res.0$Thetas, # the result from BayesALAAM
                thin=1, # should be set so that SACF is sufficiently low, 
                # important for Confidence Intervals
                tabname=NULL) # the name appended to the table that is saved

## you can improve the mixing by using a better proposal covariance
Propsigma <- cov(res.0$Thetas)
## which can be used as an argument PropSigma to BayesALAAM. This proposal 
## variance (covariance) matrix, directly regulates how big jumps we are 
## proposing, as discussed above in the section on ESS.
res.1 <- BayesALAAM(y = data$vertices$ranking50,           # dependent variable
                    ADJ = adj,           # network
                    covariates = covs,   # covariates
                    directed = TRUE,    # directed / undirected network
                    Iterations = 5000,   # number of iterations
                    saveFreq = 500,     # print and save frequency
                    PropSigma = Propsigma )

sim.0 <- get.gof.distribution(NumIterations=500, # number of vectors to draw
                              res=res.0, # the ALAAM estimation object that 
                              # contains model and results
                              burnin=100, # no. iterations discarded from 
                              # GOF distribution
                              thinning = 1000, # no. iterations between 
                              # sample points
                              contagion ='none') # should be the same as 
gof.table(obs.stats= sim.0$stats, # observed statistics included not 
          # fitted statistics
          sim.stats= sim.0$Sav.gof, # simulated goodness-of-fit statistics
          name.vec= sim.0$gof.stats.names, # names of statistics calculate, 
          # not all will be used if undirected
          tabname='ALAAMGofalt', # name of file saved
          pvalues=TRUE, # posterior predictive p-values
          # save.tab ='csv', # save a csv file or a LaTex file
          directed=TRUE)
