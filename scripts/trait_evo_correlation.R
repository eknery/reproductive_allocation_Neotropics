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

### counting pruned phylognetic trees
n_phylo = length(list.files("0_data/pruned_phylos"))

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

############################## FITTING SIMPLE MODELS #########################

### choose traits
t1 = "seed_mass"
trait = log(spp_traits[[t1]])
names(trait) = spp_traits$species

### checking normality
shapiro.test( trait )

### fit WN model
wn_fit = fitContinuous(phy= mcc_phylo, dat= trait , model = "white")
wn_aic = wn_fit$opt$aic

### ## final directory
dir_name = paste0("3_trait_results/BM/", t1)
dir_check = dir.exists(dir_name)
# create output dir if not created yet
if (dir_check == FALSE){
  dir.create(path= paste0(dir_name) )
}

### list of resutls
bm_aic_list = list()
bm_rates_list = list()

for (i in 1:n_phylo){
  ### importing phylogenetic tree
  phylo_path = paste0("0_data/pruned_phylos/pruned_phylo_", as.character(i))
  phylo = read.tree(phylo_path)
  ### fit BM model
  bm_fit = evolvcv.lite(tree = phylo,  
                    X = X 
                    )
  ### picking aicc scores and estimates
  bm_aic_list[[i]] = bm_fit$aic
  ### picking aicc scores and estimates
  bm_rates_list[[i]] = bm_fit$opt$aic
  ### check!
  print(paste0("BM fit to map: ", i))
}

### export
saveRDS(bm_aic_list, paste0(dir_name,"/bm_aic_list" ) )
saveRDS(bm_rates_list, paste0(dir_name,"/bm_rates_list" ) )

############################### FITTING R matrix ##################################

### choose traits
t1 = "seed_mass"
t2 = "sla"

### checking normality
shapiro.test( log(spp_traits[[t1]]) )
shapiro.test( log(spp_traits[[t2]]) )

### ## final directory
dir_name = paste0("3_trait_results/EVOLVCV/",t1)
dir_check = dir.exists(dir_name)
# create output dir if not created yet
if (dir_check == FALSE){
  dir.create(path= paste0("3_trait_results/EVOLVCV/",t1) )
}

### vcv matrix
X = as.matrix( log(spp_traits[,c(t1,t2)]) )
rownames(X) = spp_traits$species

### choose modelt to fit
r_models = c("1","2","3", "3b", "3c", "4")

### lists to keep model fit 
r_aic_list = list()
r_best_list = list()
r_rates_list = list()

for(i in 1:length(simmap_list) ){
  ### pick one tree
  one_simmap = simmap_list[[i]]
  ### fit vcv matrices
  r_fit = evolvcv.lite(
    tree = one_simmap$tree[[100]], 
    X = X,
    models = r_models
  )
  ### get parameter number
  k = c()
  for (m in 1:length(r_models)){
    k = c(k, r_fit[[m]]$k)
  }
  names(k) = r_models
  ### get aic scores
  aic_scores = c()
  for (m in 1:length(r_models)){
    aic_scores = c(aic_scores, r_fit[[m]]$AIC)
  }
  names(aic_scores) = r_models
  ### keep aic scores
  r_aic_list[[i]] = aic_scores
  ### delta aic
  daic = sort(unlist(aic_scores) - min(aic_scores, na.rm = T) )
  ### lowest delta aicc
  fir_model = names(daic[1])
  sec_model = names(daic[2])
  ### choosing best model
  if (daic[2] >= 2) {
    best_model = fir_model
  } 
  if (daic[2] < 2) {
    if(k[fir_model] < k[sec_model]){
      best_model = fir_model
    }
    if(k[fir_model] > k[sec_model]){
      best_model = sec_model
    }
  }
  ### keep best model and parameters
  best_model_num = which(r_models == best_model)
  r_best_list[[i]] = best_model
  r_rates_list[[i]] = r_fit[[best_model_num]]$R
  ### check
  print(paste0("VCV fit to map: ", i))
}

### export
saveRDS(r_aic_list, paste0(dir_name,"/r_aic_list.RDS") )
saveRDS(r_best_list, paste0(dir_name,"/r_best_list.RDS") )
saveRDS(r_rates_list, paste0(dir_name,"/r_rates_list.RDS") )

################################## BEST MODELS #################################

### trait name
t1 = "seed_mass"

### vcv models
dir_vcv= paste0("3_trait_results/EVOLVCV/",t1)
# aic scores
r_aic_list = readRDS(paste0(dir_vcv,"/r_aic_list.RDS") )
# best model list and parameters
r_best_list = readRDS(paste0(dir_vcv,"/r_best_list.RDS") )
r_rates_list = readRDS(paste0(dir_vcv,"/r_rates_list.RDS") )

### name of most frequent models
model_count = sort(table(unlist(r_best_list)), decreasing = T)
fir_model_name = names(model_count[1])
sec_model_name = names(model_count[2])

### indexes of most frequent models
best_model_names = unlist(r_best_list)
fir_model_index = which(best_model_names == fir_model_name)
sec_model_index = which(best_model_names == sec_model_name)

### picking correlation values
fir_cor_values = c()
for (i in fir_model_index){
  rates = r_rates_list[[i]]
  cor_value = cov2cor(rates)[1,2]
  fir_cor_values = c(fir_cor_values, cor_value)
}

### describe 
hist(fir_cor_values)
summary(fir_cor_values)
IQR(fir_cor_values)

### picking varaince values
fir_var_values = c()
for (i in fir_model_index){
  rates = r_rates_list[[i]]
  var_value = c(rates[1,1], rates[2,2])
  fir_var_values = rbind(fir_var_values, var_value)
}
fir_var_values = as.data.frame(fir_var_values)
colnames(fir_var_values) = c("seed_mass", "sla")

### describe 
summary(fir_var_values)
IQR(fir_var_values[["seed_mass"]])
IQR(fir_var_values[["sla"]])


### pick correlation values
sec_cor_values = c()
for (i in sec_model_index){
  rates = r_rates_list[[i]]
  cor_vec = c()
  for (j in 1:length(rates)){
    mtx_name = names(rates)[j]
    cor_value = cov2cor(rates[[j]])[1,2]
    cor_vec = rbind(cor_vec, c(cor_value, mtx_name) )
  }
  sec_cor_values = rbind(sec_cor_values, cor_vec)
}
### transform to dataframe
sec_cor_df = as.data.frame(sec_cor_values)
## name columns
colnames(sec_cor_df) = c("corr", "habitat_type")
### describe 
sec_cor_df %>% 
  mutate(corr = as.numeric(corr)) %>% 
  group_by(habitat_type) %>% 
  reframe(median(corr, na.rm=T), IQR(corr))

### pick correlation values
sec_var_values = c()
for (i in sec_model_index){
  rates = r_rates_list[[i]]
  var_vec = c()
  for (j in 1:length(rates)){
    mtx_name = names(rates)[j]
    var_value = c(rates[[j]][1,1] , rates[[j]][2,2] )
    var_vec = rbind(var_vec, c(var_value, mtx_name) )
  }
  sec_var_values = rbind(sec_var_values, var_vec)
}
sec_var_df = as.data.frame(sec_var_values)
## name columns
colnames(sec_var_df) = c("seed_mass", "sla", "habitat_type")
### describe 
sec_var_df %>% 
  mutate(seed_mass = as.numeric(seed_mass),
         sla = as.numeric(sla)
         ) %>% 
  group_by(habitat_type) %>% 
  reframe(median(seed_mass, na.rm=T), 
          IQR(seed_mass),
          median(sla, na.rm=T), 
          IQR(sla)
          )
