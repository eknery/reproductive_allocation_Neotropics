if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("Hmisc")) install.packages("Hmisc"); library("Hmisc")
if (!require("PupillometryR")) install.packages("PupillometryR"); library("PupillometryR")
if (!require("plyr")) install.packages("plyr"); library("plyr")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("diversitree")) install.packages("diversitree"); library("diversitree")

### create directory for pgls models
# check if dir exists
dir_check = dir.exists(paths="4_graphics")
# create dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "4_graphics")
}

### phylogenetic tree 
mcc_phylo = read.tree("0_data/pruned_mcc_phylo.nwk")

### importing habita range
habitat_range = readRDS("1_habitat_results/habitat_range.RDS")

### loading trait data
trait_mtx = read.table("0_data/trait_matrix.csv", 
                       h=T, sep=",", na.strings = "na")

############################# processing data ##############################

### n tips and nodes
n_tips = Ntip(mcc_phylo)
n_inner_nodes = mcc_phylo$Nnode

### sampled species
sampled_sp = unique(trait_mtx$species)

### habitat states
spp_states = habitat_range$range
names(spp_states) =  habitat_range$species

## state colors
state_cols = c("orange", "darkgoldenrod","darkgreen")
names(state_cols) = levels(habitat_range$range)

############################### trait analyses ##########################

### traits per species
spp_traits = trait_mtx %>% 
  group_by(species) %>% 
  reframe(height = mean(plant_height, na.rm=T) ,
          sla =  mean(sla, na.rm=T) ,
          seed = mean(seed_weight, na.rm=T) ,
          n = n()
  )

### setting regime df
species = spp_traits$species
regime = habitat_range$range

## trait name
trait_name = "sla"
## trait values
trait = spp_traits[[trait_name]]
se = spp_traits[[trait_name]] / sqrt(spp_traits[["n"]])

### species regime and traits
sp_regime_trait = data.frame(species, regime, trait, se)

### describe
sp_regime_trait %>% 
  group_by(regime) %>% 
  reframe(mean(trait), sd(trait))

### graphical param
## y axis name
if(trait_name == "height"){
  y_axis_name = "plant height (m)"
}
if(trait_name == "sla"){
  y_axis_name = "SLA (mm2/mg)"
}
if(trait_name == "seed"){
  y_axis_name = "seed mass (mg)"
}

## file name
file_name = paste0("4_graphics/",trait_name, ".tiff")

### plotting 
tiff(file_name, 
     units="in", width=2.7, height= 2, res=600)

  ggplot(data= sp_regime_trait, 
         aes(x=regime,
             y=trait, 
             color = regime, 
             fill=regime)
         ) +
    
    geom_point(position = position_jitter(width = 0.10), 
               size = 1, 
               alpha = 0.65
               ) +
    
    geom_boxplot(weight = 0.5/length(unique(sp_regime_trait$regime)),
                 linewidth = 0.2, 
                 alpha = 0.25, 
                 outlier.shape = NA
                 )+
    
    scale_fill_manual(values=state_cols)+
    scale_colour_manual(values=state_cols)+
    
    xlab("habitat range")+ ylab(y_axis_name)+
    
    theme(panel.background=element_rect(fill="white"), 
          panel.grid=element_line(colour=NULL),
          panel.border=element_rect(fill=NA,colour="black"),
          axis.title=element_text(size=8,face="bold"),
          axis.text=element_text(size=6),
          legend.position = "none")

dev.off()


############################### fitting SIMMAP #################################

### fitting ER
er_fit = fitDiscrete(phy = mcc_phylo , 
                     dat = spp_states,
                     model="ER")
### fitting SYM
sym_fit = fitDiscrete(phy = mcc_phylo , 
                      dat = spp_states,
                      model="SYM")
### fitting ARD
ard_fit = fitDiscrete(phy = mcc_phylo , 
                      dat = spp_states,
                      model="ARD")

### picking AICc scores
aicc = c(er_fit$opt$aicc, sym_fit$opt$aicc, ard_fit$opt$aicc)
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

### infer simmaps
all_maps = phytools::make.simmap(tree = mcc_phylo, 
                                 x = spp_states, 
                                 model= model,
                                 nsim=100
                                 )

### describe maps
des_map =  phytools::describe.simmap(all_maps)

### ancestral states probs
ace = des_map$ace

### setting states
# tip states probs
tip_states_probs = ace[(1+n_inner_nodes):(n_inner_nodes +n_tips), ]
tip_states_probs = tip_states_probs[,names(state_cols)]
# ancestral state probs
inner_node_probs = ace[1:n_inner_nodes,]
inner_node_probs = inner_node_probs[,names(state_cols)]

tiff("4_graphics/simmap.tiff", units="in", width=3.3, height=6, res=600)
  phytools::plotTree(tree=mcc_phylo,
                     fsize=0.6, 
                     ftype="i")
  tiplabels(pie=tip_states_probs, 
            piecol=state_cols, 
            cex=0.45)
  nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes), 
             pie= inner_node_probs,
             piecol=state_cols, 
             cex=1.2)
  axisPhylo(pos=c(0.5), font=3, cex.axis=0.5)
dev.off()

