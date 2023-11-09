if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("diversitree")) install.packages("diversitree"); library("diversitree")
if (!require("data.table")) install.packages("data.table"); require("data.table")

################################## MY FUNCTIONS #################################

source("scripts/function_calculate_aicc.R")

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
  
  bisse_full = constrain(bisse_full, q01 ~ q10)
  
  ### full bisse function w random states
  # names(rand_states) = sample(names(eco_states))
  # 
  # bisse_rand = make.bisse(tree = one_phylo, 
  #                         states = rand_states, 
  #                         sampling.f= 0.85)
  
  ### null bisse function
  bisse_null = constrain(bisse_full, lambda0 ~ lambda1, mu0 ~ mu1)
  
  ### starting
  start_bisse = starting.point.bisse(one_phylo)
  
  ### FULL MODEL find parameter values
  fit_full = find.mle(func = bisse_full,
                      x.init = start_bisse[-6]
                      )
  
  ### RANDOM MODEL find parameter values
  # fit_rand = find.mle(func = bisse_rand,
  #                     x.init = start_bisse
  #                     )
  
  ### NULL MODEL find parameter values
  fit_null = find.mle(func = bisse_null,
                      x.init = start_bisse[-c(1,3,6)]
                      )
  
  ### calculating goodness of  fit
  aicc_full = calculate_aicc(lnlik = fit_full$lnLik,
                            k = length(fit_full$par), 
                            n = n_tips
                            )
               
  # aicc_rand = calculate_aicc(lnlik = fit_rand$lnLik,
  #                            k = length(fit_rand$par),
  #                            n = n_tips
  #                            )   
  
  aicc_null = calculate_aicc(lnlik = fit_null$lnLik,
                             k = length(fit_null$par), 
                             n = n_tips
                             ) 
  
  ### vector with aicc scores             
  aicc_scores = c(aicc_full, aicc_null)
  names(aicc_scores) = c("full", "null")
  
  ### adding to list
  aicc_scores_list[[i]] = aicc_scores
  
  ### adding parameter values
  full_param_list[[i]] = fit_full$par
  #rand_param_list[[i]] = fit_rand$par
  null_param_list[[i]] = fit_null$par
  
  print(paste0("Diverfication estimates done: ", i) )
  
}

### from list to df
aicc_scores_df = data.frame( t( sapply(aicc_scores_list,c) ) )
full_param_df = data.frame( t( sapply(full_param_list,c) ) )
null_param_df = data.frame( t( sapply(null_param_list,c) ) )

### exporting AICc scores and parameters
write.table(aicc_scores_df, 
            "3_diversification_analyses/aicc_scores_df.csv",
            sep=",",
            row.names = F)

write.table(full_param_df, 
            "3_diversification_analyses/full_param_df.csv",
            sep=",",
            row.names = F)

write.table(null_param_df, 
            "3_diversification_analyses/null_param_df.csv",
            sep=",",
            row.names = F)

##################### Contrasting model and parameters ######################

### read model fit and parameter estimates
aicc_scores_df = read.table("3_diversification_analyses/q equall/aicc_scores_df.csv", 
                            sep=",", h = T)

full_param_df = read.table("3_diversification_analyses/q equall/full_param_df.csv",
                           sep=",", h = T)

null_param_df = read.table("3_diversification_analyses/q equall/null_param_df.csv",
                          sep=",", h = T)

### best fit per tree
best_models = colnames(aicc_scores_df)[apply(aicc_scores_df,1,which.min)]

### describing parameters
summary(null_param_df$lambda1 - null_param_df$mu1)

### plotting probibilities
diversi_plot = ggplot(null_param_df) +
  
  geom_density(aes(x= lambda1 - mu1),
               alpha = 0.75,
               fill ="gray",
               color = "gray",
               linetype = "solid") +
  
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
tiff("3_diversification_analyses/diversi_plot.tiff", 
     units="cm", width=7, height=6.2, res=600)
diversi_plot
dev.off()
