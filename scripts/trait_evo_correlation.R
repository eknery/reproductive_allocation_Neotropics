if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("phytools")) install.packages("phytools"); library("phytools")

### my function
source("scripts/function_calculate_aicc.R")
############################### LOADING DATA #################################

### loading phylogenetic tree
mcc_phylo = read.tree("0_data/pruned_mcc_phylo.nwk")

### loading trait data
trait_mtx = read.table("0_data/trait_matrix.csv", 
                       h=T, sep=",", na.strings = "na")

### loading occurrence count per domain
habitat_range = readRDS("1_habitat_results/habitat_range.RDS")

### loading reconstructions fro ancestral range
simmap_list = readRDS("2_reconstruction_results/simmap/simmap_list.RDS")

############################### DATA PROCESSING ###############################

### node and tip numbers
n_node = mcc_phylo$Nnode
n_tip = Ntip(mcc_phylo)

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

############################## RANGE RECONSTRUCTION ###########################

### defining states
spp_states = habitat_range$range
names(spp_states) = habitat_range$species

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


### model fit list
model_fit = list()

model_fit$er= geiger::fitDiscrete(phy = mcc_phylo , 
                                  dat = spp_states,
                                  model= models$er)

model_fit$or_sym = geiger::fitDiscrete(phy = mcc_phylo , 
                                       dat = spp_states,
                                       model= models$or_sym)

model_fit$or_asym = geiger::fitDiscrete(phy = mcc_phylo , 
                                        dat = spp_states,
                                        model= models$or_asym)

model_fit$un_sym = geiger::fitDiscrete(phy = mcc_phylo , 
                                       dat = spp_states,
                                       model= models$un_sym)

model_fit$un_asym = geiger::fitDiscrete(phy = mcc_phylo , 
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

### choosing best transition model
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
trans_param = model_fit[[best_model]]$opt

### infer simmaps
simmaps = phytools::make.simmap(tree = mcc_phylo, 
                                 x = spp_states, 
                                 model= models[[best_model]],
                                 pi = c(0.25,0.5,0.25),
                                 nsim= 100
)       

############################### TESTING CORRELATION ###########################

### choose traits
t1 = "sla"
t2 = "seed_mass"

### ## final directory
dir_name = paste0("3_trait_results/EVOLVCV/",t2)
### name of the final results
final_list_name = paste(t1,t2,"rates", sep="_")

### criando repositório para os testes OUWIE
dir_check = dir.exists(dir_name)
# create output dir if not created yet
if (dir_check == FALSE){
  dir.create(path= paste0("3_trait_results/EVOLVCV/",t2) )
}

### trait matrix
X = as.matrix(spp_traits[,c(t1,t2)])
rownames(X) = spp_traits$species

### lists to keep best model and estimates
best_model_list = list()
best_rates_list = list()

for(i in 1:length(simmaps) ){
  
  i = 1
  
  ### pick one tree
  one_simmap = simmap_list[[i]]
  ### fit vcv matrices
  vcv_fit = evolvcv.lite(
    tree = one_simmap$tree[[100]], 
    X = X
  )
  ### vector with parameter numbers
  k = c()
  ##
  k$model1 = vcv_fit$model1$k
  k$model2 = vcv_fit$model2$k
  k$model3 = vcv_fit$model3$k
  k$model4 = vcv_fit$model4$k
  k = unlist(k)
  ### vector with aicc scores
  aicc = c()
  ## model1
  aicc$model1 = calculate_aicc(lnlik = vcv_fit$model1$logLik,
                               k = vcv_fit$model1$k , 
                               n = n_tip)
  ## model2
  aicc$model2 = calculate_aicc(lnlik = vcv_fit$model2$logLik,
                               k = vcv_fit$model2$k , 
                               n = n_tip)
  ## model3
  aicc$model3 = calculate_aicc(lnlik = vcv_fit$model3$logLik,
                               k = vcv_fit$model3$k , 
                               n = n_tip)
  ## model4
  aicc$model4 = calculate_aicc(lnlik = vcv_fit$model4$logLik,
                               k = vcv_fit$model4$k , 
                               n = n_tip)
  ### delta aicc
  daicc = sort(unlist(aicc)  - min(unlist(aicc)) )
  ### lowest delta aicc
  fir_model = names(daicc[1])
  sec_model = names(daicc[2])
  ### choosing best model
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
  ### keep best model and parameters
  best_model_list[[i]] = best_model
  best_rates_list[[i]] = vcv_fit[[best_model]]$R
  ### check
  print(paste0("VCV fit to map: ", i))
}

### export
saveRDS(best_model_list, paste0(dir_name,"/best_model_list.RDS") )
saveRDS(best_rates_list, paste0(dir_name,"/best_rates_list.RDS") )

################################## BEST MODEL ##################################

### trait name
t2 = "seed_mass"

### directory name
dir_name = paste0("3_trait_results/EVOLVCV/",t2)

### best model list and parameters
best_model_list = readRDS(paste0(dir_name,"/best_model_list.RDS") )
### best model list and parameters
best_rates_list = readRDS(paste0(dir_name,"/best_rates_list.RDS") )

### pick most frequent model
model_count = table(unlist(best_model_list))
best_freq_model = names(model_count[max(model_count) == model_count])

### picking correlation values
cor_values = c()
for (i in 1:length(best_rates_list)){
  cor_mtx = cov2cor(best_rates_list[[i]])
  cor_values = c(cor_values, cor_mtx[1,2])
}

hist(cor_values)


