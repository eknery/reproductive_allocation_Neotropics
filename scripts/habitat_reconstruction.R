### libraries
if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("Hmisc")) install.packages("Hmisc"); library("Hmisc")
if (!require("PupillometryR")) install.packages("PupillometryR"); library("PupillometryR")
if (!require("plyr")) install.packages("plyr"); library("plyr")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("diversitree")) install.packages("diversitree"); library("diversitree")
if (!require("hisse")) install.packages("hisse");library("hisse")
if (!require("FAmle")) install.packages("FAmle");library("FAmle")

### loading phylogenetic tree
mcc_phylo = ape::read.tree("0_data/pruned_mcc_phylo.nwk")

### counting pruned phylognetic trees
n_phylo = length(list.files("0_data/pruned_phylos"))

### importing habitat range
habitat_range = readRDS("1_habitat_results/habitat_range.RDS")

### my function
source("scripts/function_fit_geohisse.R")

############################# PROCESSING KEY ARGUMENTS #########################

### n tips and nodes
n_tips = Ntip(mcc_phylo)
n_inner_nodes = mcc_phylo$Nnode

## "0" is the widespread area '01'
## "1" is endemic area '00' 
## "2" is endemic area '11'
### define states
species = habitat_range$species
geo_states = as.character(habitat_range$range)
geo_states[geo_states == "generalist"] = "0"
geo_states[geo_states == "forest_specialist"] = "1"
geo_states[geo_states == "open_specialist"] = "2"
geo_states = as.data.frame(cbind(species, geo_states))

### sampling fraction
sampling_f = c(0.75, 0.75, 0.75)

### states code
## explicit names
state_names = c("(00)", "(11)", "(01)")
## hidden state
histate_names = c("(00A)", "(11A)", "(01A)", "(00B)", "(11B)", "(01B)")

### state colors 
state_cols = c("lightskyblue", "mediumorchid","lightsalmon") 
names(state_cols) = levels(habitat_range$range)

############################# TRANSITION MATRICES ##############################

# trans_h1 = TransMatMakerGeoHiSSE(hidden.traits=1)

### transition matrix
trans_1 = matrix(c(NA,0,1,
                   0,NA,2,
                   0,0,NA), 
                 nrow= 3, 
                 byrow = T,
                 dimnames = list(state_names, state_names)
          )

trans_2 = matrix(c(NA,0,1,
                   0,NA,2,
                   1,2,NA), 
                 nrow= 3, 
                 byrow = T,
                 dimnames = list(state_names, state_names)
)

### transition matrix with hidden states
trans_h1 = matrix(c(NA,0,1,2,NA,NA,
                   0,NA,1,NA,2,NA,
                   0,0,NA,NA,NA,2,
                   2,NA,NA,NA,0,1,
                   NA,2,NA,0,NA,1,
                   NA,NA,2,0,0,NA), 
                 nrow= 6, 
                 byrow = T,
                 dimnames = list(histate_names, histate_names)
          )

trans_h2 = matrix(c(NA,0,1,2,NA,NA,
                    0,NA,1,NA,2,NA,
                    1,1,NA,NA,NA,2,
                    2,NA,NA,NA,0,1,
                    NA,2,NA,0,NA,1,
                    NA,NA,2,1,1,NA), 
                  nrow= 6, 
                  byrow = T,
                  dimnames = list(histate_names, histate_names)
)

############################## FITTING MODELS ##################################

### list to keep all models
model_list = list()

### model 0 - dispersal only
model_list$model_0 = GeoHiSSE(
  phy = mcc_phylo, 
  data = geo_states, 
  f= sampling_f, 
  turnover=c(1,1,0), 
  eps=c(1,1), 
  hidden.states= FALSE,
  trans.rate = trans_2, 
  root.type="madfitz",
  sann.its=1000
)

### model 1 - range dependent
model_list$model_1 = GeoHiSSE(
  phy = mcc_phylo, 
  data = geo_states, 
  f= sampling_f, 
  turnover=c(1,2,3), 
  eps=c(1,1), 
  hidden.states= FALSE,
  trans.rate = trans_2, 
  root.type="madfitz",
  sann.its=1000
)

### model 2 - hidden diversification
model_list$model_2 = GeoHiSSE(
  phy = mcc_phylo, 
  data = geo_states, 
  f= sampling_f, 
  turnover=c(1,1,0,2,3,0), 
  eps=c(1,1,1,1), 
  hidden.states= T,
  trans.rate = trans_h2, 
  root.type="madfitz",
  sann.its=1000
)

### model 3 - range dependent with hidden state
model_list$model_3 = GeoHiSSE(
  phy = mcc_phylo, 
  data = geo_states, 
  f= sampling_f, 
  turnover=c(1,2,3,4,5,6), 
  eps=c(1,1,1,1), 
  hidden.states= T,
  trans.rate = trans_h2, 
  root.type="madfitz",
  sann.its=1000
)

### AICc weights
aicc_w = GetAICWeights(list(model_0 = model_list$model_0,
                            model_1 = model_list$model_1,
                            model_2 = model_list$model_2,
                            model_3 = model_list$model_3
                            ),
                      criterion="AICc")

### name of the best model 
best_model_name = names(aicc_w[aicc_w == max(aicc_w)])
### pick the best model
best_model = model_list[[best_model_name]]

############################## MCC RECONSTRUCTION ##############################

### reconstruction
best_rec = MarginReconGeoSSE(
  phy = best_model$phy,
  data = best_model$data, 
  f = best_model$f, 
  pars = best_model$solution, 
  hidden.states= 1, 
  root.type="madfitz", 
  root.p= best_model$root.p, 
  AIC= best_model$AIC, 
  verbose=TRUE
  )

### setting states
## tip states probs
tip_states_probs = best_rec$tip.mat[,c("(00A)", "(01A)", "(11A)")]
## ancestral state probs
inner_node_probs = best_rec$node.mat[,c("(00A)", "(01A)", "(11A)")]
colnames(inner_node_probs) = levels(habitat_range$range)
## vector with ancestral states
anc_states = colnames(inner_node_probs)[apply(inner_node_probs,1,which.max)]

### plot mcc reconstruction
phytools::plotTree(tree=mcc_phylo,
                   fsize=0.6, 
                   ftype="i")
tiplabels(pie=tip_states_probs, 
          piecol=state_cols, 
          cex=0.45)
nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes), 
           pie= inner_node_probs,
           piecol=state_cols, 
           cex=1.3)
axisPhylo(pos=c(0.5), font=3, cex.axis=0.5)

############################# ALL RECONSTRUCTIONS ##############################

### vector with best model names
all_best_models = c()
### list of parameters for each phylo tree
best_param_list= list()
### list with ancestral reconstructions
anc_state_list = list()

### loop across phylo trees
for(i in 2:n_phylo){
  ### importing phylogenetic tree
  phylo_path = paste0("0_data/pruned_phylos/pruned_phylo_", as.character(i))
  phylo = read.tree(phylo_path)
  ### find best geohisse model
  best_model = fit_geohisse(phy = phylo,
               data = geo_states,
               trans_simple = trans_2,
               trans_hidden = trans_h2,
               sann= F
               )
  ### keep best model
  all_best_models = c(all_best_models, names(best_model) )
  best_param_list[[i]] = best_model[[1]]
  best_model = best_model[[1]]
  ### reconstruction
  rec = MarginReconGeoSSE(
    phy = best_model$phy,
    data = best_model$data, 
    f = best_model$f, 
    pars = best_model$solution, 
    hidden.states= 1, 
    root.type="madfitz", 
    root.p= best_model$root.p, 
    AIC= best_model$AIC, 
    verbose=TRUE
  )
  ## ancestral state probs
  inner_node_probs = rec$node.mat[,c("(00A)", "(01A)", "(11A)")]
  colnames(inner_node_probs) = levels(habitat_range$range)
  ## vector with ancestral states
  anc_states = colnames(inner_node_probs)[apply(inner_node_probs,1,which.max)]
  ### adding to list
  anc_state_list[[i]] = anc_states
  ### progress!!
  print(paste0("Ancestal reconstruction done: ", i))
}

### check results
table(all_best_models)

### exporting the name of best models
saveRDS(all_best_models, "2_reconstruction_results/all_best_models.RDS")
### exporting paramerter values
saveRDS(best_param_list, "2_reconstruction_results/best_param_list.RDS")
### exporting ancetral state list
saveRDS(anc_state_list, "2_reconstruction_results/anc_state_list.RDS")


