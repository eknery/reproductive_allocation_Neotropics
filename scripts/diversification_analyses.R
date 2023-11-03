if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("diversitree")) install.packages("diversitree"); library("diversitree")

############################### LOADING DATA #################################

### loading phylogenetic tree
mcc_phylo = read.tree("0_data/mcc_phylo.nwk")
phylo_trees = read.tree("0_data/100_rand_phylos.nwk")

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

### define states threshold
high_ths = 0.90
low_ths = (1 - high_ths)
eco_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
eco_states[af_percentage >= high_ths] = 1
eco_states[af_percentage < high_ths] = 0

## state vector
names(eco_states) = spp_count_domain$species

## vector for random states 
rand_states = eco_states 

############################### fitting SSE models ##########################

### n tips
n_tips = length(mcc_phylo$tip.label)

### list to save aicc scores
aicc_scores_list = list()

### list to save parameters
full_param_list = list()
rand_param_list = list()
null_param_list = list()

for (i in 1:length(phylo_trees) ){

  ### picking one phylogenetic tree
  one_phylo = phylo_trees[[i]]

  ### full bisse function
  bisse_full = make.bisse(tree = one_phylo, 
                          states = eco_states,
                          sampling.f= 0.85)
  
  ### full bisse function w random states
  names(rand_states) = sample(names(eco_states))
  
  bisse_rand = make.bisse(tree = one_phylo, 
                          states = rand_states, 
                          sampling.f= 0.85)
  
  ### null bisse function
  bisse_null = constrain(bisse_full, lambda0 ~ lambda1, mu0 ~ mu1)
  
  ### starting
  start_bisse = starting.point.bisse(one_phylo)
  
  ### FULL MODEL find parameter values
  fit_full = find.mle(func = bisse_full,
                      x.init = start_bisse
                      )
  
  ### RANDOM MODEL find parameter values
  fit_rand = find.mle(func = bisse_rand,
                      x.init = start_bisse
                      )
  
  ### NULL MODEL find parameter values
  fit_null = find.mle(func = bisse_null,
                      x.init = start_bisse[-c(1,3)]
                      )
  
  ### calculating goodness of  fit
  aicc_full = calculate_aicc(lnlik = fit_full$lnLik,
                            k = length(fit_full$par), 
                            n = n_tips
                            )
               
  aicc_rand = calculate_aicc(lnlik = fit_rand$lnLik,
                             k = length(fit_rand$par),
                             n = n_tips
                             )   
  
  aicc_null = calculate_aicc(lnlik = fit_null$lnLik,
                             k = length(fit_null$par), 
                             n = n_tips
                             ) 
  
  ### vector with aicc scores             
  aicc_scores = c(aicc_full, aicc_rand, aicc_null)
  names(aicc_scores) = c("full", "rand", "null")
  
  ### adding to list
  aicc_scores_list[[i]] = aicc_scores
  
  ### adding parameter values
  full_param_list[[i]] = fit_full$par
  rand_param_list[[i]] = fit_rand$par
  null_param_list[[i]] = fit_null$par
  
  print(paste0("Diverfication estimates done: ", i) )
  
}

##################### Contrasting model and parameters ######################

### read model and parameter estimates
mcmc_full = read.table("3_diversification_analyses/mcmc_full.csv", sep=",", h = T)
mcmc_rand = read.table("3_diversification_analyses/mcmc_rand.csv", sep=",", h = T)
mcmc_null = read.table("3_diversification_analyses/mcmc_null.csv", sep=",", h = T)

### describing probabilities and parameters
median(mcmc_null$p)
IQR(mcmc_null$p)


median(mcmc_full$mu1)
IQR(mcmc_full$mu1)

### plotting probibilities
prob_plot = ggplot() +
  
  geom_density(data = mcmc_null,
               aes(x= p), 
               alpha = 0.75,
               fill ="white", 
               color = "black",
               linetype = "dashed") +
  
  geom_density(data = mcmc_full,
               aes(x= p),
               alpha = 0.75,
               fill ="black", 
               color = "black") +
  
  scale_x_continuous(limits = c(-184,-162), 
                     expand = c(0, 0)) +
  
  xlab("Posterior probability") +
  ylab("Density")  +
  
  theme(panel.background= element_rect(fill="white"),
        panel.grid= element_line(colour=NULL),
        panel.border= element_rect(fill=NA,colour="black"),
        axis.title= element_text(size= 10, face="bold"),
        axis.text.x= element_text(size= 8, angle = 0),
        axis.text.y= element_text(size= 8, angle = 0),
        legend.position= "none")

## exporting
tiff("3_diversification_analyses/prob_plot.tiff", 
     units="cm", width=7, height=6.2, res=600)
prob_plot
dev.off()


### plotting lambda
lambda_plot = ggplot(data = mcmc_full) +
  
  geom_density(aes(x= lambda0), 
               alpha = 0.75,
               fill ="orange", 
               color = "orange") +
  geom_density(aes(x= lambda1),
               alpha = 0.75,
               fill ="darkgreen", 
               color = "darkgreen") +
  
  scale_x_continuous(limits = c(0,1.1), 
                     expand = c(0, 0)) +
  
  xlab("Lambda") +
  ylab("Density")  +
  
  theme(panel.background= element_rect(fill="white"),
        panel.grid= element_line(colour=NULL),
        panel.border= element_rect(fill=NA,colour="black"),
        axis.title= element_text(size= 10, face="bold"),
        axis.text.x= element_text(size= 8, angle = 0),
        axis.text.y= element_text(size= 8, angle = 0),
        legend.position= "none")

## exporting
tiff("3_diversification_analyses/lambda_plot.tiff", 
     units="cm", width=7, height=6.2, res=600)
lambda_plot
dev.off()

### plotting MU
mu_plot = ggplot(data = mcmc_full) +
  
  geom_density(aes(x= mu0), 
               alpha = 0.75,
               fill ="orange", 
               color = "orange") +
  geom_density(aes(x= mu1),
               alpha = 0.75,
               fill ="darkgreen", 
               color = "darkgreen") +
  
  scale_x_continuous(limits = c(-0.25,0.75), 
                     expand = c(0, 0)) +
  
  xlab("Mu") +
  ylab("Density")  +
  
  theme(panel.background= element_rect(fill="white"),
        panel.grid= element_line(colour=NULL),
        panel.border= element_rect(fill=NA,colour="black"),
        axis.title= element_text(size= 10, face="bold"),
        axis.text.x= element_text(size= 8, angle = 0),
        axis.text.y= element_text(size= 8, angle = 0),
        legend.position= "none")

## exporting
tiff("3_diversification_analyses/mu_plot.tiff", 
     units="cm", width=7, height=6.2, res=600)
mu_plot
dev.off()


### plotting diversification rates
diversification_plot = ggplot(data = mcmc_full) +
  
  geom_density(aes(x= lambda0 - mu0), 
               alpha = 0.75,
               fill ="orange", 
               color = "orange") +
  geom_density(aes(x= lambda1 - mu1),
               alpha = 0.75,
               fill ="darkgreen", 
               color = "darkgreen") +
  
  scale_x_continuous(limits = c(-0.3,0.9), 
                     expand = c(0, 0)) +
  
  xlab("Diversification rates") +
  ylab("Density")  +
  
  theme(panel.background= element_rect(fill="white"),
        panel.grid= element_line(colour=NULL),
        panel.border= element_rect(fill=NA,colour="black"),
        axis.title= element_text(size= 10, face="bold"),
        axis.text.x= element_text(size= 8, angle = 0),
        axis.text.y= element_text(size= 8, angle = 0),
        legend.position= "none")

## exporting
tiff("3_diversification_analyses/diversification_plot.tiff", 
     units="cm", width=7, height=6.2, res=600)
diversification_plot
dev.off()


### plotting transition rates
q_plot = ggplot(data = mcmc_full) +
  
  geom_density(aes(x= q01), 
               alpha = 0.75,
               fill ="orange", 
               color = "orange") +
  geom_density(aes(x= q10),
               alpha = 0.75,
               fill ="darkgreen", 
               color = "darkgreen") +
  
  scale_x_continuous(limits = c(-0.1,0.7), 
                     expand = c(0, 0)) +
  
  xlab("Transition rates") +
  ylab("Density")  +
  
  theme(panel.background= element_rect(fill="white"),
        panel.grid= element_line(colour=NULL),
        panel.border= element_rect(fill=NA,colour="black"),
        axis.title= element_text(size= 10, face="bold"),
        axis.text.x= element_text(size= 8, angle = 0),
        axis.text.y= element_text(size= 8, angle = 0),
        legend.position= "none")

## exporting
tiff("3_diversification_analyses/q_plot.tiff", 
     units="cm", width=7, height=6.2, res=600)
q_plot
dev.off()
