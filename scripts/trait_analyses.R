if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("OUwie")) install.packages("OUwie"); library("OUwie")

################################ MY FUNCTIONS #################################

source("scripts/function_fit_evo_models.R")
source("scripts/function_choose_best.R")

############################### LOADING DATA #################################

### RANDOM SEED
set.seed(42)

### loading phylogenetic tree
pruned_mcc_phylo = read.tree("0_data/pruned_mcc_phylo.nwk")

### counting pruned phylognetic trees
n_phylo = length(list.files("0_data/pruned_phylos"))

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

### loading trait data
trait_mtx = read.table("0_data/trait_matrix.csv", 
                       h=T, sep=",", na.strings = "na")

############################### DATA PROCESSING ###############################

### sampled species
sampled_sp = unique(trait_mtx$species)

### define ecological state
high_ths = 0.90
low_ths = (1 - high_ths)
eco_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
eco_states[af_percentage >= high_ths] = "specialist"
eco_states[af_percentage < high_ths] = "generalist"

names(eco_states) = spp_count_domain$species

### keeping only sampled species
keep_eco_states = eco_states[names(eco_states) %in% sampled_sp]

############################### EDA ########################################

### trait values per species
sp_traits = trait_mtx %>% 
  group_by(species) %>% 
  reframe(height = mean(plant_height, na.rm=T) ,
          sla =  mean(sla, na.rm=T) ,
          seed = mean(seed_weight, na.rm=T) ,
          n = n()
  )

############################ ANCESTRAL STATE RECONSTRUCTION ##############################

### ancestral node numbers
n_node = pruned_mcc_phylo$Nnode

### list of ancestral states for each phylo
q_values_list = list()

### list of ancestral states for each phylo
anc_states_list = list()

for (i in 1:n_phylo){
  
  ### importing phylogenetic tree
  pruned_phylo_path = paste0("0_data/pruned_phylos/pruned_phylo_", as.character(i))
  pruned_phylo = read.tree(pruned_phylo_path)
    
  ### fitting equal rates
  er_fit = fitDiscrete(phy = pruned_phylo , 
                        dat = keep_eco_states,
                        model="ER")
  
  ### fitting symmetric
  ard_fit = fitDiscrete(phy = pruned_phylo , 
                        dat =keep_eco_states,
                        model="ARD")
  
  ### picking AICc scores
  aicc= c(er_fit$opt$aicc, ard_fit$opt$aicc)
  names(aicc) = c("ER","ARD")
  
  ### chossing best transition model
  if (aicc[["ER"]] <= aicc[["ARD"]]) {
    model = "ER"
  } else {
    model = "ARD"
  }
  
  ### choosing Q matrix
  if(model == "ER"){
    q_values= c("ER", er_fit$opt$q12, er_fit$opt$q21)
    names(q_values) = c("model","q12","q21")
  }
  if(model == "ARD"){
    q_values= c("ARD",ard_fit$opt$q12, ard_fit$opt$q21)
    names(q_values) = c("model", "q12","q21")
  }
  
  ### keeping best fit Q matrix
  q_values_list[[i]] = q_values
  
  ### prior for ancestral state
  pi = c(1,0)
  names(pi) = c("generalist", "specialist")
  
  ### infer simmaps
  all_maps = phytools::make.simmap(tree = pruned_phylo, 
                                   x = keep_eco_states, 
                                   model= model,
                                   nsim=100,
                                   pi = pi
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
  
  print(paste0("Ancestal reconstruction done:", i))
  
} 

### exporting ancetral state list
saveRDS(q_values_list, "2_trait_analyses/q_values_list.RDS")

### exporting ancetral state list
saveRDS(anc_states_list, "2_trait_analyses/anc_states_list.RDS")

########################## comparing transition models ######################

### importing Q values
q_values_list = readRDS("2_trait_analyses/q_values_list.RDS")

### from list to df
q_values_df = data.frame( t( sapply(q_values_list,c) ) )

### most common best fit
table(q_values_df$model)

### transition rates
q12 = as.numeric(q_values_df$q12)
q21 = as.numeric(q_values_df$q21)

### describing transition rates
mean(q12)
sd(q12)

mean(q21)
sd(q21)
################################ OUWIE ########################################

### loading ancetral state list
anc_states_list = readRDS("2_trait_analyses/anc_states_list.RDS")

### setting regime df
species = sp_traits$species
regime = keep_eco_states

## trait name
trait_name = "sla"
## trait values
trait = sp_traits[[trait_name]]
se = sp_traits[[trait_name]] / sqrt(sp_traits[["n"]])

sp_regime_trait = data.frame(species, regime, trait, se)

dir_check = dir.exists(paths=paste("2_trait_analyses/OUWIE/",trait_name, sep="") )
# create output dir if not created yet
if (dir_check == FALSE){
  dir.create(path= paste("2_trait_analyses/OUWIE/",trait_name, sep="") )
}

### model to fit 
all_models = c("BM1","BMS","OU1", "OUM","OUMV")

### model fit table
all_best_models = data.frame(matrix(NA, nrow= n_phylo, ncol=4))
colnames(all_best_models) = c(c("model","llik","aicc","delta_aicc"))

### best estiamtes list
all_best_estimates = list()

for (i in 1:n_phylo){
  
  ### importing phylogenetic tree
  pruned_phylo_path = paste0("0_data/pruned_phylos/pruned_phylo_", as.character(i))
  pruned_phylo = read.tree(pruned_phylo_path)
  
  ### adding ancestral states to phylo tree
  pruned_phylo$node.label = anc_states_list[[i]]
  
  
  all_fits = fit_evo_models(phy= pruned_phylo, 
                            data= sp_regime_trait,
                            mserr = 'known',
                            models_to_fit = all_models)
  
  best_choice = choose_best(all_fits)
  
  all_best_models[i,] = best_choice$best_fit
  all_best_estimates[[i]] = best_choice$best_estimates
  
  print(paste0("Trait evolution done: ", i) )
  
}

### export path
exp_path = paste0("2_trait_analyses/OUWIE/", trait_name)

### exporting model fit
saveRDS(all_best_models, paste0(exp_path, "/all_best_models.RDS") )

### exporting best estimates list
saveRDS(all_best_estimates, paste0(exp_path, "/all_best_estimates.RDS") )

######################## Comparing models and parameters #######################

### choose a trait!
trait_name = "height"

### import path
imp_path = paste0("2_trait_analyses/OUWIE/", trait_name)

### importing model fit
all_best_models = readRDS( paste0(imp_path, "/all_best_models.RDS"))

### importing model parameters
all_best_estimates =  readRDS( paste0(imp_path, "/all_best_estimates.RDS"))

### most frequent best fit 
count_models = table(all_best_models$model)
most_freq_model = names(count_models)[count_models == max(count_models)]

### vector w position of best fit
i_best = which(all_best_models$model == most_freq_model)

### select parameter estimates from most freq best model
freq_best_estimates = all_best_estimates[i_best]

### from list to df
est_df = data.frame( t( sapply(freq_best_estimates,c) ) )


apply(est_df, MARGIN = 2, mean)
apply(est_df, MARGIN = 2, sd)
