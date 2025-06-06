if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("phytools")) install.packages("phytools"); library("phytools")
if (!require("nlme")) install.packages("nlme"); library("nlme")
if (!require("OUwie")) install.packages("OUwie"); library("OUwie")

################################ MY FUNCTIONS #################################

source("scripts/function_fit_evo_models.R")
source("scripts/function_choose_best.R")

############################### LOADING DATA #################################

### loading phylogenetic tree
mcc_phylo = read.tree("0_data/pruned_mcc_phylo.nwk")

### counting pruned phylognetic trees
n_phylo = length(list.files("0_data/pruned_phylos"))

### loading trait data
trait_mtx = read.table("0_data/trait_matrix.csv", 
                       h=T, sep=",", na.strings = "na")

### loading occurrence count per domain
habitat_range = readRDS("1_habitat_results/habitat_range.RDS")

### loading list of ancestral states
anc_state_list = readRDS("2_reconstruction_results/geohisse/anc_state_list.RDS")

############################### DATA PROCESSING ###############################

### sampled species
sampled_sp = unique(trait_mtx$species)

### defininf states
spp_states = habitat_range$range
names(spp_states) = habitat_range$species

### trait values per species
spp_traits = trait_mtx %>% 
  mutate(
    seed_wei_mg = fruit_weight_mg/seed_number,
  ) %>% 
  group_by(species) %>% 
  reframe(
    sla =  median(leaf_sla, na.rm=T),
    seed_mass = median(seed_wei_mg),
    n = n()
  )

##################################### OUWIE ####################################

### criando repositório para os testes OUWIE
dir_check = dir.exists(paths="3_trait_results/OUWIE" )
# create output dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "3_trait_results/OUWIE" )
}

### setting regime df
species = spp_traits$species
regime = spp_states

## trait name
trait_name = "rel_inflor"
## trait values
trait =  spp_traits[[trait_name]] 
se = sd(trait) / sqrt(spp_traits[["n"]])

## ouwie table
sp_regime_trait = data.frame(species, regime, trait, se)

dir_check = dir.exists(paths=paste("3_trait_results/OUWIE/",trait_name, sep="") )
# create output dir if not created yet
if (dir_check == FALSE){
  dir.create(path= paste0("3_trait_results/OUWIE/",trait_name) )
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
  phylo_path = paste0("0_data/pruned_phylos/pruned_phylo_", as.character(i))
  phylo = read.tree(phylo_path)
  
  ### adding ancestral states to phylo tree
  phylo$node.label = anc_state_list[[i]]
  
  ### fitting all models
  all_fits = fit_evo_models(phy= phylo, 
                            data= sp_regime_trait,
                            #mserr = 'known',
                            models_to_fit = all_models)
  ### chose best model
  best_choice = choose_best(all_fits)
  ### take best estimates
  all_best_models[i,] = best_choice$best_fit
  all_best_estimates[[i]] = best_choice$best_estimates
  
  print(paste0("Trait evolution, tree done: ", i) )
  
}

### export path
exp_path = paste0("3_trait_results/OUWIE/", trait_name)
### exporting model fit
saveRDS(all_best_models, paste0(exp_path, "/all_best_models.RDS") )
### exporting best estimates list
saveRDS(all_best_estimates, paste0(exp_path, "/all_best_estimates.RDS") )
 
