copurchase <- read.csv("copurchase.csv")
products <- read.csv("products.csv")

# Question 1
#filter products to only books that have a salesrank of <= 150000 and salesrank !=-1
new_product <- subset(products, products$group == 'Book' & products$salesrank <= 150 & products$salesrank != -1)
book_id <- unique(new_product$id)
new_copurchase <- subset(copurchase, copurchase$Source %in% book_id & copurchase$Target %in% book_id)


library(statnet)


# View the first rows of the edgelist to make sure it imported correctly:
head(new_copurchase)
# Convert the edgelist to a network object in statnet format:
advice <- as.network.matrix(new_copurchase, matrix.type = "edgelist") 
advice

set.vertex.attribute(advice, "rating",new_product$rating)
set.vertex.attribute(advice, "salesrank",new_product$salesrank)
set.vertex.attribute(advice, "review_cnt",new_product$review_cnt)
set.vertex.attribute(advice, "downloads",new_product$downloads)



options(ergm.loglik.warn_dyads=FALSE) #Whether or not a warning should be issued when sample space constraints render the observed number of dyads ill-defined

# Ergm Terms are statistics: They are some deterministic function of the ties, node attributes, and edge covariates of a network.
help("ergm-terms",package = "ergm") # Documentation that contains definitions for all of the terms we are using
# ex. what does "mutual" test and how is it calculated
# We will use the ergm-terms to perform hypothesis testing using ERGMs
# But we can note that any of the ERGM terms can also be examined directly for your observed network, by creating a formula in R

# Look at Endogenous statistics: terms based on only ties in the advice network
summary(advice ~ edges)                     # Number of edges (ties)
summary(advice ~ mutual) 
summary(advice ~ odegree(0:5))              # Outdegree distribution. (# of nodes with outdegree of 0, # nodes outdegree of 1, etc.)
# Remember, respondents could nominate at most five employees in our survey
summary(advice ~ idegree(0:65))             # Indegree distribution.
summary(advice ~ gwodegree(log(2),fixed=T)) # One parameter summarizing outdegree distribution - tendency against outdegree hubs
summary(advice ~ gwidegree(log(2),fixed=T)) # One parameters summarizing indegree distribution - tendency against indegree hubs
summary(advice ~ desp(1:5))                 # Pairs of nodes with one shared partner, two shared partners, etc.
summary(advice ~ dgwesp(log(2),fixed = T))


model1 <- ergm(advice ~ edges                 # This is  a tendency towards a greater number of advice ties existing. Based on a statistic counting the number of ties.
               # Structural patterns
               + mutual                      # This is a tendency towards reciprocity for the advice ties. Based on a statistic counting the number of reciprocated ties.
                  # This is the effect of every 100 messages sent from i->j on likelihood of an advice tie. Based on a weighted sum of advice ties x 100s of messages sent
               
               # Model constraints
               , constraints =~ bd(maxout=5) # This constraint enforces the maximum outdegree is 5
) 
summary(model1) 
summary(advice ~ nodecov("downloads"))
summary(advice ~ nodecov("salesrank"))
summary(advice ~ nodecov("review_cnt"))
summary(advice ~ nodecov("rating"))

model2 <- ergm(advice ~ 
               #edges                 # This is  a tendency towards a greater number of advice ties existing. Based on a statistic counting the number of ties.
               # Structural patterns
               + mutual                      # This is a tendency towards reciprocity for the advice ties. Based on a statistic counting the number of reciprocated ties.
               # This is the effect of every 100 messages sent from i->j on likelihood of an advice tie. Based on a weighted sum of advice ties x 100s of messages sent
               + gwidegree(log(2), fixed = T)                 # Inverted preferential attachment (indegree)
               + gwodegree(2, fixed = T, cutoff = 5)              # Inverted preferential attachment (outdegree)
               + dgwesp(log(2), type = "OTP", fixed = T, cutoff =5)
               + nodecov('rating')
               + nodecov('salesrank')
               + nodecov('review_cnt')
               + nodecov('downloads')
               # Model constraints
               , constraints =~ bd(maxout=5)
               , control = control.ergm(MCMC.effectiveSize = 50)# This constraint enforces the maximum outdegree is 5
) 
summary(model2) 

library(texreg)
screenreg(list("model1"=model1,"model2"=model2))
# Screenshot this output and include it in your report

# Ergm uses a random algorithm, so every time you estimate the model you may get slightly different results
# Thus, you may want to save your R environment, into a .Rdata file:
# save.image("Lab2_files.RData")
# And then you can re-load your results into R later by clicking to open the file, or by running the line of code below:
# load("Lab2_files.RData")

# -------------------------------------------------------------------------------------------------
######################## PART III: Model diagnostics ########################
# -------------------------------------------------------------------------------------------------
# MCMC diagnostics - will save to a pdf in the current working directory

pdf('model1diagnostics.pdf')              # Open a pdf file to save to
mcmc.diagnostics(model1) # Run the markov chain monte carlo diagnostics
dev.off()                # Closes and saves the pdf

pdf('model2diagnostics.pdf')              # Open a pdf file to save to
mcmc.diagnostics(model2) # Run the markov chain monte carlo diagnostics
dev.off()                # Closes and saves the pdf

# -------------------------------------------------------------------------------------------------
# Goodness of fit test - will display in RStudio
# Check how well the estimated model captures certain features of the observed network, for example triangles in the network.
# -------------------------------------------------------------------------------------------------
# Look at networks simulated according to model 1
# This first command simulates 100 networks.
# These networks, if we use sufficient burnin steps in the markov chain used to generate them,
# may be thought of as random samples from the joint probability distribution that is our fitted ERGM.
library('igraph')
net_layout <- layout_with_fr(advice_igraph)
sim1 <- simulate(model1, burnin=100000, interval=100000, nsim=100, verbose=T)  # Uses the ergm model to simulate a null model
# Plot the first of the simulated networks
sim1_net1 <- igraph::graph.adjacency(as.matrix.network(sim1[[1]]))
igraph::plot.igraph(sim1_net1,layout=net_layout,edge.color="brown",  
                    vertex.color = 'grey',edge.arrow.size=.4)                                                               
# Plot the 10th simulated network
sim1_net10 <- igraph::graph.adjacency(as.matrix.network(sim1[[10]]))
igraph::plot.igraph(sim1_net10,layout=net_layout,edge.color="red",  
                    vertex.color = 'grey',edge.arrow.size=.4)

# Repeat, now looking at networks simulated according to model 2
sim2 <- simulate(model2, burnin=100000, interval=100000, nsim=100, verbose=T)  # Uses the ergm model to simulate a null model
# Plot the first of the simulated networks
sim2_net1 <- igraph::graph.adjacency(as.matrix.network(sim2[[1]]))
igraph::plot.igraph(sim2_net1,layout=net_layout,edge.color="grey",  
                    vertex.color = 'grey',edge.arrow.size=.4)                                                               
# Plot the 10th simulated network
sim2_net10 <- igraph::graph.adjacency(as.matrix.network(sim2[[10]]))
igraph::plot.igraph(sim2_net10,layout=net_layout,edge.color="purple",  
                    vertex.color = 'grey',edge.arrow.size=.4)

# -------------------------------------------------------------------------------------------------
# Extract the number of triangles from each of the 100 samples and
# compare the distribution of triangles in the sampled networks with the observed network
# -------------------------------------------------------------------------------------------------
# Model 1:
model1.tridist <- sapply(1:100, function(x) summary(sim1[[x]] ~triangle)) # Extracts the triangle data from the simulated networks
hist(model1.tridist,xlim=c(0,1000),breaks=10)                             # Plots that triangle distribution as a histogram, change xlim to change the x-axis range if necessary
advice.tri <- summary(advice ~ triangle)                                    # Stores the number of observed triangles
advice.tri
arrows(advice.tri,20, advice.tri, 0.5, col="red", lwd=3)                      # Adds an arrow to the plotted histogram
c(obs=advice.tri,mean=mean(model1.tridist),sd=sd(model1.tridist),
  tstat=abs(mean(model1.tridist)-advice.tri)/sd(model1.tridist))

# Model 2:
model2.tridist <- sapply(1:100, function(x) summary(sim2[[x]] ~triangle)) # Extracts the triangle data from the simulated networks
hist(model2.tridist,xlim=c(0,1000),breaks=10)                             # Plots that triangle distribution as a histogram, change xlim to change the x-axis range if necessary
arrows(advice.tri,20, advice.tri, 0.5, col="red", lwd=3)                    # Adds an arrow to the plotted histogram
c(obs=advice.tri,mean=mean(model2.tridist),sd=sd(model2.tridist),
  tstat=abs(mean(model2.tridist)-advice.tri)/sd(model2.tridist))

# -------------------------------------------------------------------------------------------------
# Test the goodness of fit of the model
# Compiles statistics for these simulations as well as the observed network, and calculates p-values 
# -------------------------------------------------------------------------------------------------

# Model 1:
# It may take a second for this command to run.
gof1 <- gof(model1, verbose=T, burnin=1e+5, interval=1e+5, control = control.gof.ergm(nsim = 200))
# If you run below and then wouldn't see the plot, trypar(mar=c(2,2,2,2))
dev.off()           # Clear any other plots from the plot window
plot(gof1)          # Plot the goodness of fit
# Note: This should produce five separate plots that you should look through.
#       In RStudio, scroll between the plots using the arrow buttons
gof1                # Display the goodness of fit info in the console

# Model 2:
# It may take a second for this command to run.
gof2 <- gof(model2, verbose=T, burnin=1e+5, interval=1e+5, control = control.gof.ergm(nsim = 200))
dev.off()
plot(gof2)
gof2
