### libraries
if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("Hmisc")) install.packages("Hmisc"); library("Hmisc")
if (!require("PupillometryR")) install.packages("PupillometryR"); library("PupillometryR")
if (!require("plyr")) install.packages("plyr"); library("plyr")

### loading phylogenetic tree
mcc_phylo = read.tree("0_data/pruned_mcc_phylo.nwk")

### counting pruned phylognetic trees
n_phylo = length(list.files("0_data/pruned_phylos"))

### read data
habitat_mtx = read.table("0_data/habitat_matrix.csv", 
                         h=T, sep=",", na.strings = "na")

##################################### EDA ####################################

apply(habitat_mtx[,-1], FUN = sum, MARGIN = 2)

################################ CLASSIFICATION ##############################

### habitat names
habitat_names = habitat_mtx %>% 
  select(! species) %>% 
  colnames()

### habitats
forest = habitat_names[grepl(pattern= "Floresta", habitat_names)]
non_forest = habitat_names[!grepl(pattern= "Floresta", habitat_names)]

### presence per habitat
habitat_pres = habitat_mtx %>% 
  pivot_longer(cols = any_of(habitat_names), 
               names_to = "habitat",
               values_to = "presence")

### assigning habitat types
habitat_type = habitat_pres %>% 
  mutate(type = case_when(
    (habitat %in% forest)     & presence == 1  ~ "forest",
    (habitat %in% non_forest) & presence == 1  ~ "non_forest"
    )
  ) %>% 
  filter( !is.na(type) )

### counting habitat types
habitat_count = habitat_type %>% 
  pivot_wider(names_from = type, 
              values_from = presence,
              values_fill = 0) %>% 
  group_by(species) %>% 
  reframe(
          n_forest = sum(forest), 
          n_open = sum(non_forest)
          )

### defining habtiat range
habitat_range = habitat_count %>% 
  mutate(
    range = case_when(
      n_forest == 0  & n_open >= 1 ~ "open-specialist",
      n_forest >= 1  & n_open >= 1 ~ "open-generalist",
      n_forest >= 1  & n_open == 0 ~ "forest-generalist"
    )
  ) %>% 
  mutate(range = factor(range, levels = c("open-specialist",
                                          "open-generalist",
                                          "forest-generalist")
                        )
         )

table(habitat_range$range)

### with state probibility
habitat_range_prob = habitat_range %>% 
    mutate(prob = 1) %>% 
    pivot_wider(names_from = range,
                names_expand = T,
                values_from = prob,
                values_fill = 0
                )
  
### exporting 
saveRDS(habitat_range, "1_habitat_results/habitat_range.RDS")

######################## ANCESTRAL STATE RECONSTRUCTION ########################

### importing habita range
habitat_range = readRDS("1_habitat_results/habitat_range.RDS")

### defininf states
spp_states = habitat_range$range
names(spp_states) = habitat_range$species

### ancestral node numbers
n_node = mcc_phylo$Nnode

### list of ancestral states for each phylo
q_values_list = list()

### list of ancestral states for each phylo
anc_states_list = list()

for (i in 1:n_phylo){
  
  ### importing phylogenetic tree
  phylo_path = paste0("0_data/pruned_phylos/pruned_phylo_", as.character(i))
  phylo = read.tree(phylo_path)
  
  ### fitting ER
  er_fit = fitDiscrete(phy = phylo , 
                       dat = spp_states,
                       model="ER")
  ### fitting SYM
  sym_fit = fitDiscrete(phy = phylo , 
                        dat = spp_states,
                        model="SYM")
  ### fitting ARD
  ard_fit = fitDiscrete(phy = phylo , 
                        dat = spp_states,
                        model="ARD")
  
  ### picking AICc scores
  aicc= c(er_fit$opt$aicc, sym_fit$opt$aicc, ard_fit$opt$aicc)
  names(aicc) = c("ER","SYM","ARD")
  
  ### parameter numbers
  k = c(1, 3, 6)
  names(k) = c("ER","SYM","ARD")
  
  ### delta aicc
  daicc = sort(aicc  - min(aicc))
  
  ### lowest daicc
  fir_model = names(daicc[1])
  sec_model = names(daicc[2])
  
  ### chossing best transition model
  if (fir_model == "ER") {
    model = "ER"
  } 
  if (fir_model != "ER"){
    if(daicc[sec_model] >= 2){
      model = fir_model
    } 
    if(daicc[2] < 2 & k[fir_model] < k[sec_model]){
      model = fir_model
    }
    if(daicc[2] < 2 & k[fir_model] > k[sec_model]){
      model = sec_model
    }
  }
  
  ### choosing Q matrix
  if(model == "ER"){
    q_values = unlist(er_fit$opt[ which(grepl("q", names(er_fit$opt))) ] )
    names(model) = "model"
    q_values = c(model, q_values)
  }
  if(model == "SYM"){
    q_values = unlist(sym_fit$opt[ which(grepl("q", names(sym_fit$opt))) ] )
    names(model) = "model"
    q_values = c(model, q_values)
  }
  if(model == "ARD"){
    q_values = unlist(ard_fit$opt[ which(grepl("q", names(ard_fit$opt))) ] )
    names(model) = "model"
    q_values = c(model, q_values)
  }
  
  ### keeping best fit Q matrix
  q_values_list[[i]] = q_values
  
  ### infer simmaps
  all_maps = phytools::make.simmap(tree = phylo, 
                                   x = spp_states, 
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
  
  print(paste0("Ancestal reconstruction done: ", i))
  
} 

### exporting q-values list
saveRDS(q_values_list, "1_habitat_results/q_values_list.RDS")

### exporting ancetral state list
saveRDS(anc_states_list, "1_habitat_results/anc_states_list.RDS")

########################## comparing transition models ######################

### importing Q values
q_values_list = readRDS("1_habitat_results/q_values_list.RDS")

### from list to df
q_values_df = data.frame( t( sapply(q_values_list,c) ) )

### most common best fit
tab_models = table(q_values_df$model)
best_model= names(tab_models[tab_models == max(tab_models)])

### transition rates
q_values_df %>% 
  filter(model == best_model) %>% 
  select(!model) %>% 
  as.numeric() %>% 
  summary()
