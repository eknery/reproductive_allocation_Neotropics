if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("phytools")) install.packages("phytools"); library("phytools")
if (!require("slouch")) install.packages("slouch"); library("slouch")


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
    rel_inflor = 100*inflorescence_length_cm/(plant_height_m*100)
  ) %>% 
  group_by(species) %>% 
  reframe(
    plant_hei = median(plant_height_m, na.rm=T) ,
    sla =  median(leaf_sla, na.rm=T),
    inflor_len = median(inflorescence_length_cm, na.rm=T),
    rel_inflor = median(rel_inflor, na.rm=T),
    fruit_wei = median(fruit_weight_mg, na.rm=T) ,
    seed_num = median(seed_number, na.rm=T) ,
    seed_wei = median(seed_wei_mg),
    n = n()
  )

################################### SLOUCH #####################################

### reordering
spp_traits = spp_traits[match(mcc_phylo$tip.label, spp_traits$species), ]
## check
spp_traits$species == mcc_phylo$tip.label

### pick a trait
trait = "seed_wei"

### fit slouch
slouch_model = slouch.fit(
  phy = mcc_phylo,
  species = spp_traits$species,
  response = log(spp_traits[[trait]] )
)

slouch_model$V
