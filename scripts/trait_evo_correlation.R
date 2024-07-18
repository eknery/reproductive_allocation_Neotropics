if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("phytools")) install.packages("phytools"); library("phytools")

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
  mutate(seed_wei_mg = fruit_weight_mg/seed_number) %>% 
  group_by(species) %>% 
  reframe(
    plant_hei = median(plant_height_m, na.rm=T) ,
    leaf_sla =  median(leaf_sla, na.rm=T) ,
    inflor_len = median(inflorescence_length_cm, na.rm=T),
    fruit_wei = median(fruit_weight_mg, na.rm=T) ,
    seed_num = median(seed_number, na.rm=T) ,
    seed_wei = median(seed_wei_mg),
    n = n()
  )

evolvcv.lite(sunfish.tree,sunfish.data[,2:3])