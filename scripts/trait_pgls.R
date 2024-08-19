if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("phytools")) install.packages("phytools"); library("phytools")
if (!require("nlme")) install.packages("nlme"); library("nlme")

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

################################## PGLS #########################################

### choosing a trait
trait_name = "seed_mass"
trait = log(spp_traits[[trait_name]])
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
