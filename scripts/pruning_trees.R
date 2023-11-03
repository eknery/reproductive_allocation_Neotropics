### require packages
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("phytools")) install.packages("phytools"); library("phytools")
if (!require("geiger")) install.packages("geiger"); library("geiger")

### load mcc phylogenetic tree
mcc_phylo = read.tree("0_data/mcc_phylo.nwk")

### loading random sample of trees
rand_phylos = read.tree("0_data/100_rand_phylos.nwk")

### loading trait data
trait_mtx = read.table("0_data/trait_matrix.csv", 
                       h=T, sep=",", na.strings = "na")
### sampled species
sampled_species = unique(trait_mtx$species)

##################################### PRUNING #############################

### prunning mcc phylogenetic tree
pruned_mcc_phylo = drop.tip(mcc_phylo, mcc_phylo$tip.label[-match(sampled_species, mcc_phylo$tip.label)])

### exporting prunned mcc
write.tree(pruned_mcc_phylo, file = "0_data/pruned_mcc_phylo.nwk", append = FALSE, digits = 10, tree.names = FALSE)

dir_check = dir.exists(paths="0_data/pruned_phylos" )
# create output dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "0_data/pruned_phylos" )
}

### pruning trees to sampled species
pruned_phylos = rand_phylos
for (i in 1:length(rand_phylos)){
  pruned_phylos[[i]] = drop.tip(rand_phylos[[i]], 
                                rand_phylos[[i]]$tip.label[-match(sampled_species, rand_phylos[[i]]$tip.label)])
}

### exporting phylogenetic trees
for (i in 1:length(pruned_phylos)){
  write.tree(pruned_phylos[[i]], 
             file = paste("0_data/pruned_phylos/pruned_phylo_", as.character(i), sep=""), 
             append = FALSE, digits = 10, tree.names = FALSE)
}
