### libraries
if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("Hmisc")) install.packages("Hmisc"); library("Hmisc")
if (!require("PupillometryR")) install.packages("PupillometryR"); library("PupillometryR")
if (!require("plyr")) install.packages("plyr"); library("plyr")
if (!require("ape")) install.packages("ape"); library("ape")

### loading phylogenetic tree
mcc_phylo = ape::read.tree("0_data/pruned_mcc_phylo.nwk")

### counting pruned phylognetic trees
n_phylo = length(list.files("0_data/pruned_phylos"))

### importing habitat range
habitat_range = readRDS("1_habitat_results/habitat_range.RDS")

############################## Ancestral reconstruction ########################

### defining states
spp_states = habitat_range$range
names(spp_states) = habitat_range$species

### ancestral node numbers
n_node = mcc_phylo$Nnode

### my models
models = list()

models$er = matrix(c(0,1,1,
                     1,0,1,
                     1,1,0), 3, byrow = T)


models$or_sym =  matrix(c(0,1,0,
                          1,0,2,
                          0,2,0), 3, byrow = T)

models$or_asym  = matrix(c(0,1,0,
                           2,0,3,
                           0,4,0), 3, byrow = T)

models$un_sym = matrix(c(0,1,2,
                         1,0,3,
                         2,3,0), 3, byrow = T)

models$un_asym = matrix(c(0,1,2,
                          3,0,4,
                          5,6,0), 3, byrow = T)
### parameter numbers
k = c("er" = 1,
      "or_sym" = 2, 
      "or_asym" = 4,
      "un_sym" = 3,
      "un_asym" = 6
)
### list of best parameters for each phylo
param_list = list()
### list of simmaps for each phylo
simmap_list = list()
### list of ancestral states for each phylo
anc_state_list = list()

for (i in 1:n_phylo){
  
  ### importing phylogenetic tree
  phylo_path = paste0("0_data/pruned_phylos/pruned_phylo_", as.character(i))
  phylo = read.tree(phylo_path)
  
  ### model fit list
  model_fit = list()
 
  model_fit$er= geiger::fitDiscrete(phy = phylo , 
                       dat = spp_states,
                       model= models$er)
  
  model_fit$or_sym = geiger::fitDiscrete(phy = phylo , 
                                  dat = spp_states,
                                  model= models$or_sym)
  
  model_fit$or_asym = geiger::fitDiscrete(phy = phylo , 
                                   dat = spp_states,
                                   model= models$or_asym)
  
  model_fit$un_sym = geiger::fitDiscrete(phy = phylo , 
                                     dat = spp_states,
                                     model= models$un_sym)
  
  model_fit$un_asym = geiger::fitDiscrete(phy = phylo , 
                                dat = spp_states,
                                model= models$un_asym)
  
  ### picking AICc scores
  aicc = c("er" = model_fit$er$opt$aicc, 
           "or_sym" = model_fit$or_sym$opt$aicc,
           "or_asym" = model_fit$or_asym$opt$aicc,
           "un_sym" = model_fit$un_sym$opt$aicc,
           "un_asym" = model_fit$un_asym$opt$aicc
  )
  
  ### delta aicc
  daicc = sort(aicc  - min(aicc))
  
  ### lowest daicc
  fir_model = names(daicc[1])
  sec_model = names(daicc[2])
  
  ### chossing best transition model
  if (daicc[2] >= 2) {
    best_model = fir_model
  } 
  if (daicc[2] < 2) {
    if(k[fir_model] < k[sec_model]){
      best_model = fir_model
    }
    if(k[fir_model] > k[sec_model]){
      best_model = sec_model
    }
  }
  
  ### keep best model parameters
  param_list[[i]] = rbind(model_fit[[best_model]]$opt)
  names(param_list)[i] = best_model
  
  ### infer simmaps
  all_maps = phytools::make.simmap(tree = phylo, 
                                   x = spp_states, 
                                   model= models[[best_model]],
                                   pi = c(0.25,0.5,0.25),
                                   nsim= 100
                                   )                               
  
  ### describe maps
  des_map =  phytools::describe.simmap(all_maps)
  ### keep map
  simmap_list[[i]] = des_map
  ### ancestral states probs
  ace = des_map$ace
  ### all states
  all_states = colnames(ace)[apply(ace,1,which.max)]
  ### ancestral node states
  anc_states = all_states[1:n_node]
  ### adding to list
  anc_state_list[[i]] = anc_states
  
  print(paste0("Ancestal reconstruction done: ", i))
  
} 

### exporting q-values list
saveRDS(param_list, "2_reconstruction_results/simmap/param_list.RDS")
### exporting q-values list
saveRDS(simmap_list, "2_reconstruction_results/simmap/simmap_list.RDS")
### exporting ancetral state list
saveRDS(anc_state_list, "2_reconstruction_results/simmap/anc_states_list.RDS")

########################## comparing transition models ######################

### importing Q values
param_list = readRDS("2_reconstruction_results/simmap/param_list.RDS")

### from list to df
param_df = data.frame( t( sapply(param_list,c) ) )
colnames(param_df) = colnames(param_list[[1]])
row.names(param_df) = NULL
param_df$model = colnames(sapply(param_list,c))

### most common best fit
tab_models = table(param_df$model)
best_model = names(tab_models[tab_models == max(tab_models)])

### pick one paramenter
param_name = "q12"
param_vector = unlist(param_df[param_name])

### outlier boundaries
med = median(param_vector, na.rm = T) 
bound= IQR(param_vector, na.rm = T)*1.5
up_bound = med + bound
dw_bound = med - bound
  
### cleaning data
clean_param_vector = param_vector[param_vector < up_bound & param_vector > dw_bound]

### describe
median(clean_param_vector); IQR(clean_param_vector)

hist(clean_param_vector)
