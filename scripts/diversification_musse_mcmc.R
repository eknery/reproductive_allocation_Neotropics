################################# packages ##################################

if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("diversitree")) install.packages("diversitree"); library("diversitree")
if (!require("data.table")) install.packages("data.table"); require("data.table")

############################### LOADING DATA #################################

### loading phylogenetic tree
mcc_phylo = read.tree("0_data/pruned_mcc_phylo.nwk")

### loading occurrence count per domain
habitat_range = readRDS("1_habitat_results/habitat_range_alt.RDS")

########################## processing data ############################

### N nodes
n_nodes = mcc_phylo$Nnode

### habitat states
spp_states = habitat_range$range
names(spp_states) =  habitat_range$species

### states to integers
spp_int = spp_states
spp_int[spp_states == "open-specialist"] = 1
spp_int[spp_states == "open-generalist"] = 2
spp_int[spp_states == "outside"] = 3
spp_int = as.integer(spp_int)
names(spp_int) = names(spp_states)

###  number of states
k = length(unique(spp_int))

################################# MUSSE ##############################

### full bisse function
musse_full = make.musse(tree = mcc_phylo, 
                        states = spp_int,
                        k = k,
                        sampling.f= 0.75)

### starting points
start_musse = starting.point.musse(mcc_phylo, k = k)

### bounds
up_bounds = (start_musse + 0.1)  * (n_nodes) 
lw_bounds = rep(0, length(start_musse))

### mcmcm search
mcmc_full = mcmc(lik = musse_full,
                 x.init = start_musse,
                 lower = lw_bounds,
                 upper = up_bounds,
                 w = 0.1,
                 nsteps = 100
                 )


### null bisse function
musse_null = constrain(musse_full, 
                       lambda2 ~ lambda1, 
                       lambda3 ~ lambda1,
                       mu2 ~ mu1, 
                       mu3 ~ mu1)

null_params = c("lambda1", "mu1", "q12", "q13", "q21", "q23", "q31", "q32")

### start point and bounds
start_null = start_musse[names(start_musse) %in% null_params]
up_null = up_bounds[names(up_bounds) %in% null_params]
lw_null = lw_bounds[names(up_bounds) %in% null_params]

mcmc_null = mcmc(lik = musse_null,
                 x.init = start_null,
                 upper = up_null,
                 lower = lw_null,
                 w = 0.1,
                 nsteps = 100)


ggplot() +
  geom_histogram(data = mcmc_null, 
                 aes(x= p), 
                 fill = "gray70",
                 alpha = 0.5)+
  geom_histogram(data = mcmc_full, 
                 aes(x= p), 
                 fill = "black",
                 alpha= 0.5) 

### exporting 
saveRDS(mcmc_full, "3_diversification_results/mcmc_full.RDS")
saveRDS(mcmc_null, "3_diversification_results/mcmc_null.RDS")
