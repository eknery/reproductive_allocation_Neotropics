if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("phytools")) install.packages("phytools"); library("phytools")
if (!require("OUwie")) install.packages("OUwie"); library("OUwie")
if (!require("nlme")) install.packages("nlme"); library("nlme")

################################ MY FUNCTIONS #################################

source("scripts/function_fit_evo_models.R")
source("scripts/function_choose_best.R")

############################### LOADING DATA #################################

### loading phylogenetic tree
mcc_phylo = read.tree("0_data/pruned_mcc_phylo.nwk")

### counting pruned phylognetic trees
n_phylo = length(list.files("0_data/pruned_phylos"))

### loading occurrence count per domain
habitat_range = readRDS("1_habitat_results/habitat_range.RDS")

### loading trait data
trait_mtx = read.table("0_data/trait_matrix.csv", 
                       h=T, sep=",", na.strings = "na")

############################### DATA PROCESSING ###############################

### sampled species
sampled_sp = unique(trait_mtx$species)

### defininf states
spp_states = habitat_range$range
names(spp_states) = habitat_range$species

### trait values per species
spp_traits = trait_mtx %>% 
  group_by(species) %>% 
  reframe(height = mean(plant_height, na.rm=T) ,
          sla =  mean(sla, na.rm=T) ,
          seed = mean(seed_weight, na.rm=T) ,
          n = n()
  )

################################## PGLS #########################################

### choosing a trait
trait = log(spp_traits$seed)
names(trait) = spp_traits$species

### fitting models
fit_bm = fitContinuous(phy= mcc_phylo, dat = trait,  model="BM")
fit_ou = fitContinuous(phy= mcc_phylo, dat = trait,  model="OU")

### choosing model aicc
if(fit_bm$opt$aicc < fit_ou$opt$aicc){
  sigsq = fit_bm$opt$sigsq
  cor_str = corBrownian(sigsq, phy = mcc_phylo, form= ~1)
}
if(fit_bm$opt$aicc > fit_ou$opt$aicc & (fit_bm$opt$aicc - fit_ou$opt$aicc) >= 2 ){
  alpha = fit_ou$opt$alpha
  cor_str = corMartins(alpha, phy = mcc_phylo, form= ~1)
}

### fitting pgls
fit_gls = gls(trait ~ spp_states,
              correlation= cor_str, 
              method = "REML")

summary(fit_gls)
plot(fit_gls)

### checking residuals
res = resid(fit_gls)[1:nrow(spp_traits)]
hist(res)
shapiro.test(res)

##################################### OUWIE ####################################

### criando repositório para os testes OUWIE
dir_check = dir.exists(paths="2_trait_results/OUWIE" )
# create output dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "2_trait_results/OUWIE" )
}

### loading ancetral state list
anc_states_list = readRDS("1_habitat_results/anc_states_list.RDS")

### setting regime df
species = spp_traits$species
regime = spp_states

## trait name
trait_name = "seed"
## trait values
trait = log( spp_traits[[trait_name]] )
se = sd(trait) / sqrt(spp_traits[["n"]])

## ouwie table
sp_regime_trait = data.frame(species, regime, trait, se)

dir_check = dir.exists(paths=paste("2_trait_results/OUWIE/",trait_name, sep="") )
# create output dir if not created yet
if (dir_check == FALSE){
  dir.create(path= paste0("2_trait_results/OUWIE/",trait_name) )
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
  phylo$node.label = anc_states_list[[i]]
  
  ### fitting all models
  all_fits = fit_evo_models(phy= phylo, 
                            data= sp_regime_trait,
                            mserr = 'known',
                            models_to_fit = all_models)
  
  best_choice = choose_best(all_fits)
  
  all_best_models[i,] = best_choice$best_fit
  all_best_estimates[[i]] = best_choice$best_estimates
  
  print(paste0("Trait evolution, tree done: ", i) )
  
}

### export path
exp_path = paste0("2_trait_results/OUWIE/", trait_name)

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
