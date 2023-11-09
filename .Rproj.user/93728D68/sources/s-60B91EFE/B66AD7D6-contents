if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("Hmisc")) install.packages("Hmisc"); library("Hmisc")
if (!require("PupillometryR")) install.packages("PupillometryR"); library("PupillometryR")
if (!require("plyr")) install.packages("plyr"); library("plyr")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")


### create directory for pgls models
# check if dir exists
dir_check = dir.exists(paths="4_graphics")
# create dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "4_graphics", showWarnings = , recursive = FALSE, mode = "0777")
}


### phylogenetic tree 
mcc_phylo = read.tree("0_data/mcc_phylo.nwk")

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

### loading trait data
trait_mtx = read.table("0_data/trait_matrix.csv", 
                       h=T, sep=",", na.strings = "na")


############################# common parameters ##############################

### n tips and nodes
n_tips = Ntip(mcc_phylo)
n_inner_nodes = mcc_phylo$Nnode

### define ecological state
high_ths = 0.90
low_ths = (1 - high_ths)
eco_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
eco_states[af_percentage >= high_ths] = "specialist"
eco_states[af_percentage < high_ths] = "generalist"

names(eco_states) = spp_count_domain$species

############################### trait analyses ##########################

### sampled species
sampled_sp = unique(trait_mtx$species)

### traits per species
sp_traits = trait_mtx %>% 
  group_by(species) %>% 
  reframe(height = mean(plant_height, na.rm=T) ,
          sla =  mean(sla, na.rm=T) ,
          seed = mean(seed_weight, na.rm=T) ,
          n = n()
  )

### keeping only sampled species
keep_eco_states = eco_states[names(eco_states) %in% sampled_sp]

### setting regime df
species = sp_traits$species
regime = keep_eco_states

## trait name
trait_name = "sla"
## trait values
trait = sp_traits[[trait_name]]
se = sp_traits[[trait_name]] / sqrt(sp_traits[["n"]])

### species regime and traits
sp_regime_trait = data.frame(species, regime, trait, se)

### changing to factor
sp_regime_trait = sp_regime_trait %>% 
  mutate(regime = factor(regime, levels= c( "specialist", "generalist")))

### describe
sp_regime_trait %>% 
  group_by(regime) %>% 
  reframe(mean(trait), sd(trait))

### graphical param

## state colors
state_cols=c( "darkorange","darkgreen")
names(state_cols)=c("generalist",  "specialist")

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
     units="in", width=2.3, height= 2, res=600)

ggplot(data= sp_regime_trait, 
       aes(x=regime,
           y=trait, 
           color = regime, 
           fill=regime)
       ) +
  
  geom_point(position = position_jitter(width = 0.07), 
             size = 2, 
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
  
  scale_x_discrete(labels=c("specialist" = "forest-specialist", 
                            "generalist" = "habitat-generalist"))+
  
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title=element_text(size=8,face="bold"),
        axis.text=element_text(size=6),
        legend.position = "none")

dev.off()


############################### fitting SIMMAP #################################


### fitting equal rates
er_fit = fitDiscrete(phy = mcc_phylo , 
                     dat =  eco_states,
                     model="ER")


### fitting symmetric
ard_fit = fitDiscrete(phy = mcc_phylo , 
                      dat = eco_states,
                      model="ARD")

### picking AICc scores
aicc= c(er_fit$opt$aicc, ard_fit$opt$aicc)
names(aicc) = c("ER","ARD")

### chossing best transition model
if (aicc[["ER"]] <= aicc[["ARD"]]) {
  model = "ER"
} else {
  model = "ARD"
}

### prior for ancestral state
pi = c(1,0)
names(pi) = c("generalist", "specialist")

### infer simmaps
all_maps = phytools::make.simmap(tree = mcc_phylo, 
                                 x = eco_states, 
                                 model= model,
                                 nsim=100,
                                 pi= pi
                                 )

### describe maps
des_map =  phytools::describe.simmap(all_maps)


### ancestral states probs
ace = des_map$ace

### setting states
# tip states probs
tip_states_probs = ace[(1+n_inner_nodes):(n_inner_nodes +n_tips), ]
# ancestral state probs
inner_node_probs = ace[1:n_inner_nodes,]
# state colors
state_cols=c( "darkorange","darkgreen")
names(state_cols)=c("generalist",  "specialist")


tiff("4_graphics/simmap.tiff", units="in", width=4, height=6, res=600)
  phytools::plotTree(tree=mcc_phylo,
                     fsize=0.6, 
                     ftype="i")
  tiplabels(pie=tip_states_probs, 
            piecol=state_cols, 
            cex=0.4)
  nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes), 
             pie= inner_node_probs,
             piecol=state_cols, 
             cex=1.1)
  axisPhylo(pos=c(0.5), font=3, cex.axis=0.5)
dev.off()


