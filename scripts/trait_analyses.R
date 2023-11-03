if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("OUwie")) install.packages("OUwie"); library("OUwie")

if (!require("nlme")) install.packages("nlme"); library("nlme") ## tirar se não for usar PGLS

################################ MY FUNCTIONS #################################

source("scripts/function_fit_evo_models.R")
source("scripts/function_choose_best.R")

############################### LOADING DATA #################################

### RANDOM SEED
set.seed(42)

### loading phylogenetic tree
pruned_mcc_phylo = read.tree("0_data/pruned_mcc_phylo.nwk")

### loading trait data
trait_mtx = read.table("0_data/trait_matrix.csv", 
                       h=T, sep=",", na.strings = "na")

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

### counting pruned phylognetic trees
n_phylo = length(list.files("0_data/pruned_phylos"))

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

### describing sampling
sp_sampling = trait_mtx %>% 
  group_by(species) %>%
  reframe(n =n())  

summary(sp_sampling$n)

### median values per species
sp_traits = trait_mtx %>% 
  group_by(species) %>% 
  reframe(height = median(plant_height, na.rm=T) ,
          sla =  median(sla, na.rm=T) ,
          seed = median(seed_weight, na.rm=T) 
  )

############################ ANCESTRAL STATE RECONSTRUCTION ##############################

### ancestral node numbers
n_node = pruned_mcc_phylo$Nnode

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
  sym_fit = fitDiscrete(phy = pruned_phylo , 
                        dat =keep_eco_states,
                        model="SYM")
  
  ### picking AICc scores
  aicc= c(er_fit$opt$aicc, sym_fit$opt$aicc)
  names(aicc) = c("ER","SYM")
  
  ### chossing best transition model
  if (aicc[["ER"]] <= aicc[["SYM"]]) {
    model = "ER"
  } else {
    model = "SYM"
  }
  
  ### infer simmaps
  all_maps = phytools::make.simmap(tree = pruned_phylo, 
                                   x = keep_eco_states, 
                                   model= model,
                                   nsim=100
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
saveRDS(anc_states_list, "2_trait_analyses/anc_states_list.RDS")

################################ OUWIE ########################################

### loading ancetral state list
anc_states_list = readRDS("2_trait_analyses/anc_states_list.RDS")

### setting regime df
species = sp_traits$species
regime = keep_eco_states

## trait name
trait_name = "seed"
## trait values
trait = sp_traits[[trait_name]]

sp_regime_trait = data.frame(species, regime, trait)

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
  
  
  all_fits = fit_evo_models(tree= pruned_phylo, 
                            regimes= sp_regime_trait, 
                            models_to_fit = all_models)
  
  best_choice = choose_best(all_fits)
  
  all_best_models[i,] = best_choice$best_fit
  all_best_estimates[[i]] = best_choice$best_estimates
  
  print(paste0("Trait evolution done: ", i) )
  
}

### export path
exp_path = paste0("2_trait_analyses/OUWIE/", trait_name)

### expor model fit
saveRDS(all_best_models, paste0(exp_path, "/all_best_models.RDS") )

### exporting best estimates list
saveRDS(all_best_estimates, paste0(exp_path, "/all_best_estimates.RDS") )
