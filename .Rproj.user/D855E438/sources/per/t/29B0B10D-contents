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


models$dir_sym =  matrix(c(0,1,0,
                           1,0,2,
                           0,2,0), 3, byrow = T)

models$dir_asym  = matrix(c(0,1,0,
                            2,0,3,
                            0,4,0), 3, byrow = T)

models$nondir_sym = matrix(c(0,1,2,
                             1,0,3,
                             2,3,0), 3, byrow = T)

models$nondir_asym = matrix(c(0,1,2,
                              3,0,4,
                              5,6,0), 3, byrow = T)
### parameter numbers
k = c("er" = 1,
      "dir_sym" = 2, 
      "dir_asym" = 4,
      "nondir_sym" = 3,
      "nondir_asym" = 6
)

### list of best parameters for each phylo
param_list = list()

### list of ancestral states for each phylo
anc_states_list = list()

for (i in 1:n_phylo){
  
  ### importing phylogenetic tree
  phylo_path = paste0("0_data/pruned_phylos/pruned_phylo_", as.character(i))
  phylo = read.tree(phylo_path)
  
  ### model fit list
  model_fit = list()
 
  model_fit$er= geiger::fitDiscrete(phy = phylo , 
                       dat = spp_states,
                       model= models$er)
  
  model_fit$dir_sym = geiger::fitDiscrete(phy = phylo , 
                                  dat = spp_states,
                                  model= models$dir_sym)
  
  model_fit$dir_asym = geiger::fitDiscrete(phy = phylo , 
                                   dat = spp_states,
                                   model= models$dir_asym)
  
  model_fit$nondir_sym = geiger::fitDiscrete(phy = phylo , 
                                     dat = spp_states,
                                     model= models$nondir_sym)
  
  model_fit$nondir_asym = geiger::fitDiscrete(phy = phylo , 
                                dat = spp_states,
                                model= models$nondir_asym)
  
  ### picking AICc scores
  aicc = c("er" = model_fit$er$opt$aicc, 
           "dir_sym" = model_fit$dir_sym$opt$aicc,
           "dir_asym" = model_fit$dir_asym$opt$aicc,
           "nondir_sym" = model_fit$nondir_sym$opt$aicc,
           "nondir_asym" = model_fit$nondir_asym$opt$aicc
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
  ### ancestral states probs
  ace = des_map$ace
  ### all states
  all_states = colnames(ace)[apply(ace,1,which.max)]
  ### ancestral node states
  anc_states = all_states[1:n_node]
  ### adding to list
  anc_states_list[[i]] = anc_states
  
  print(paste0("Ancestal reconstruction done: ", i))
  
} 

### exporting q-values list
saveRDS(param_list, "1_habitat_results/param_list.RDS")
### exporting ancetral state list
saveRDS(anc_states_list, "1_habitat_results/anc_states_list.RDS")

########################## comparing transition models ######################

### importing Q values
param_list = readRDS("1_habitat_results/param_list.RDS")

### from list to df
param_df = data.frame( t( sapply(param_list,c) ) )
colnames(param_df) = colnames(param_list[[1]])
row.names(param_df) = NULL
param_df$model = colnames(sapply(param_list,c))

### most common best fit
tab_models = table(param_df$model)
best_model= names(tab_models[tab_models == max(tab_models)])

### pick one paramenter
param_name = "q23"
param_vector = unlist(param_df[param_name])

### outlier boundaries
med = median(param_vector, na.rm = T) 
bound= IQR(param_vector, na.rm = T)*1.5
up_bound = med + bound
dw_bound = med - bound
  
### cleaning data
clean_param_vector = param_vector[param_vector < up_bound & param_vector > dw_bound]

### describe
mean(clean_param_vector); sd(clean_param_vector)
