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

############################### FITTING CORRELATION ###########################

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

### trait matrix
X = as.matrix( log(spp_traits[,c(t1,t2)]) )
rownames(X) = spp_traits$species

### choose modelt to fit
vcv_models = c("1","2","3", "3b", "3c", "4")

### lists to keep model fit 
aic_list = list()
best_model_list = list()
best_rates_list = list()

for(i in 1:length(simmaps) ){
  ### pick one tree
  one_simmap = simmap_list[[i]]
  ### fit vcv matrices
  vcv_fit = evolvcv.lite(
    tree = one_simmap$tree[[100]], 
    X = X,
    models = vcv_models
  )
  ### get parameter number
  k = c()
  for (m in 1:length(vcv_models)){
    k = c(k, vcv_fit[[m]]$k)
  }
  names(k) = vcv_models
  ### get aic scores
  aic_scores = c()
  for (m in 1:length(vcv_models)){
    aic_scores = c(aic_scores, vcv_fit[[m]]$AIC)
  }
  names(aic_scores) = vcv_models
  ### keep aic scores
  aic_list[[i]] = aic_scores
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
  best_model_num = which(vcv_models == best_model)
  best_model_list[[i]] = best_model
  best_rates_list[[i]] = vcv_fit[[best_model_num]]$R
  ### check
  print(paste0("VCV fit to map: ", i))
}

### export
saveRDS(aic_list, paste0(dir_name,"/aic_list.RDS") )
saveRDS(best_model_list, paste0(dir_name,"/best_model_list.RDS") )
saveRDS(best_rates_list, paste0(dir_name,"/best_rates_list.RDS") )

################################## BEST MODELS #################################

### trait name
t1 = "seed_mass"

### directory name
dir_name = paste0("3_trait_results/EVOLVCV/",t1)

### best model list and parameters
best_model_list = readRDS(paste0(dir_name,"/best_model_list.RDS") )
### best model list and parameters
best_rates_list = readRDS(paste0(dir_name,"/best_rates_list.RDS") )

### name of most frequent models
model_count = sort(table(unlist(best_model_list)), decreasing = T)
fir_model_name = names(model_count[1])
sec_model_name = names(model_count[2])

### indexes of most frequent models
best_model_names = unlist(best_model_list)
fir_model_index = which(best_model_names == fir_model_name)
sec_model_index = which(best_model_names == sec_model_name)

### picking correlation values
fir_cor_values = c()
for (i in fir_model_index){
  rates = best_rates_list[[i]]
  cor_value = cov2cor(rates)[1,2]
  fir_cor_values = c(fir_cor_values, cor_value)
}

### describe 
hist(fir_cor_values)
summary(fir_cor_values)

### pick correlation values
sec_cor_values = list()
loop = 1
for (i in sec_model_index){
  rates = best_rates_list[[i]]
  cor_vec = c()
  for (j in 1:length(rates)){
    cor_value = cov2cor(rates[[j]])[1,2]
    cor_vec = c(cor_vec, cor_value)
  }
  names(cor_vec) = names(best_rates_list[[i]])
  sec_cor_values[[loop]] = cor_vec
  loop = loop + 1
}

### transform to dataframe
sec_cor_df = data.frame(matrix(unlist(sec_cor_values), 
                  nrow=length(sec_cor_values), 
                  byrow=TRUE
                  )
           
           )
## name columns
colnames(sec_cor_df) = names(best_rates_list[[i]])

### export 
saveRDS(fir_cor_values, paste0(dir_name,"/fir_model_cor.RDS") )
saveRDS(sec_cor_df, paste0(dir_name,"/sec_model_cor.RDS") )
