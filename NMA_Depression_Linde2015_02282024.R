
library(meta)
library(netmeta)
data(Linde2015)
trts<-c("TCA","SSRI","SNRI","NRI", "Low-dose SARI",
        "NaSSa","rMAO-A","Hypericum","Placebo")
outcomes<-c("Earlyresponse", "Earlyremission")



#(1)Earlyresponse #
p1<-pairwise(treat=list(treatment1,treatment2,treatment3), 
             event=list(resp1,resp2,resp3),
             n=list(n1,n2,n3), 
             studlab=id,
             data=Linde2015,sm="OR") 
net1<-netmeta(p1,common=TRUE, seq=trts,ref ="Placebo",small="undesirable")
summary(net1, ref="Placebo")

library(dplyr)
data <- p1 %>% 
  filter(treat1 =="TCA" | treat1 =="Hypericum" | treat1 =="NaSSa" ) %>%
  filter(treat2 =="TCA" | treat2 =="Hypericum" | treat2 =="NaSSa" ) 
net2<-netmeta(data,common=TRUE, seq=trts<-c("TCA","Hypericum","NaSSa"),
              ref ="TCA",small="undesirable")
summary(net2, ref="TCA")

# 2b. Comparison of all treatment pairs
pairs2 <- netleague(net2, digits = 2,seq=trts)
View(pairs2$common)
forest(net2, ref="TCA", xlab="OR", pooled="common")

# network graph
netgraph(net2, 
         plastic = FALSE,
         number.of.studies = TRUE,
         multiarm = TRUE,
         seq = "optimal",
         cex = 1,
         rotate = 0,
         cex.points = n.trts,
         labels = paste0(gsub("+", " +\n", trts, 
                              fixed = TRUE),
                         "\n(n=", round(n.trts), ")"))

# Treatment Ranking
# P-scores allow ranking the treatments on a continuous 0-1 scale
# Based on frequentist point estimates and standard errors
netrank(net1, small.values = "undesirable", method = "P-score")
netrank(net1, small.values = "undesirable", method = "SUCRA")

# 2a. Forest Plot
forest(net1, ref="Placebo", xlab="OR", pooled="common")

# 2b. Comparison of all treatment pairs
pairs <- netleague(net1, digits = 2,seq=trts)
View(pairs$common)
View(pairs$random)

# 2c.Inconsistency
# global approach
decomp.design(net1)
# local approach
print(netsplit(net1), digits=2)

# Inconsistency refers to disagreement between direct and indirect evidence for a particular treatment comparison. 
# Inconsistency occurs 
# when the results of a direct comparison (i.e., a head-to-head comparison of two treatments) do not agree with
#      the results of an indirect comparison (i.e., a comparison that is inferred through a common comparator).

############################################################################################################ Bayesian NMA
############################################################################################################ Bayesian NMA
############################################################################################################ Bayesian NMA
# install JAGS 4.3.0 for R version 4.1.2
# install.packages("https://cran.r-project.org/src/contrib/Archive/rjags/rjags_4-10.tar.gz", type="Source", repos=NULL)
# Self-learning material:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6635665/
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/bayesnma.html

library(rjags)
library(gemtc)

# network setup: data with long arm-based format
NMA_Baye <- mtc.network(data.ab = data.long.arm.base)
summary(NMA_Baye)
plot(NMA_Baye)
# network model - fixed effect model
model_Baye <- mtc.model(NMA_Baye, linearModel= "fixed", n.chain=4)

# need to install JAGS software
mcmc_Baye <- mtc.run(model_Baye, n.adapt=5000, n.iter=20000, thin=1)
summary(mcmc_Baye)
plot(mcmc_Baye)
gelman.plot(mcmc_Baye)

# Rank probabilities indicate the probability for each treatment to be best, second best, etc.
ranks <- rank.probability(mcmc_Baye)
print(sucra(ranks))

# 2a. Forest Plot
forest(relative.effect(mcmc_Baye, t1="Sham"), digits=4)

# 2b. Comparison of all treatment pairs
round(relative.effect.table(mcmc_Baye, covariate=NA), 2)

# 2c.Inconsistency
# Split nodes
nodesplit <- mtc.nodesplit(NMA_Baye, linearModel = "fixed", n.adapt = 5000, n.iter = 20000, thin = 1)

# detach(package:gemtc, unload = TRUE)
