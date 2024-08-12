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

### loading trait data
trait_mtx = read.table("0_data/trait_matrix.csv", 
                       h=T, sep=",", na.strings = "na")
### importing habita range
habitat_range = readRDS("1_habitat_results/habitat_range.RDS")

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
state_cols = c("forestgreen", "brown","darkorange") 
names(state_cols) = levels(habitat_range$range)

############################### trait analyses ##########################

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

### setting regime df
species = spp_traits$species
regime = habitat_range$range

## trait name
trait_name = "seed_mass"
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
if(trait_name == "plant_hei"){
  y_axis_name = "plant height (m)"
}
if(trait_name == "sla"){
  y_axis_name = "SLA (mm2/mg)"
}
if(trait_name == "inflor_len"){
  y_axis_name = "inflorescence length (cm)"
}
if(trait_name == "seed_mass"){
  y_axis_name = "seed mass (mg)"
}

### trait plot
trait_plot  = ggplot(data= sp_regime_trait, 
         aes(x=regime,
             y= trait, 
             color = regime, 
             fill=regime)
         ) +
    
    geom_point(position = position_jitter(width = 0.10), 
               size = 1, 
               alpha = 0.75
               ) +
    
    geom_boxplot(weight = 0.5/length(unique(sp_regime_trait$regime)),
                 linewidth = 0.3, 
                 color= "black",
                 alpha = 0.50, 
                 outlier.shape = NA
                 )+
    
    scale_fill_manual(values=state_cols)+
    scale_colour_manual(values=state_cols)+
    scale_x_discrete(labels= c(
      "rainforest",
      "generalist",
      "open-vegetation"
      )
    ) +
    
    xlab("habitat type")+ ylab(y_axis_name)+
    
    theme(panel.background=element_rect(fill="white"), 
          panel.grid=element_line(colour=NULL),
          panel.border=element_rect(fill=NA,colour="black"),
          axis.title=element_text(size=8,face="bold"),
          axis.text=element_text(size=6),
          legend.position = "none")

## file name
file_name = paste0("4_graphics/",trait_name, ".tiff")
### plotting 
tiff(file_name, 
     units="cm", width=7, height= 6.5, res=600)
 print(trait_plot)
dev.off()


############################### fitting SIMMAP #################################

### my models
my_models = list()

my_models$er = matrix(c(0,1,1,
                        1,0,1,
                        1,1,0), 3, byrow = T)

my_models$dir_sym =  matrix(c(0,1,0,
                              1,0,2,
                              0,2,0), 3, byrow = T)

my_models$dir_asym  = matrix(c(0,1,0,
                               2,0,3,
                               0,4,0), 3, byrow = T)

my_models$nondir_sym = matrix(c(0,1,2,
                                1,0,3,
                                2,3,0), 3, byrow = T)

my_models$nondir_asym = matrix(c(0,1,2,
                                 3,0,4,
                                 5,6,0), 3, byrow = T)


### fit models
er_fit = fitDiscrete(phy = mcc_phylo , 
                     dat = spp_states,
                     model= my_models$er)

dir_sym_fit = fitDiscrete(phy = mcc_phylo , 
                          dat = spp_states,
                          model= my_models$dir_sym)

dir_asym_fit = fitDiscrete(phy = mcc_phylo , 
                           dat = spp_states,
                           model= my_models$dir_asym)

nondir_sym_fit = fitDiscrete(phy = mcc_phylo , 
                      dat = spp_states,
                      model= my_models$nondir_sym)

nondir_asym_fit = fitDiscrete(phy = mcc_phylo , 
                      dat = spp_states,
                      model= my_models$nondir_asym)

### picking AICc scores
aicc = c("er" = er_fit$opt$aicc, 
         "nondir_sym" = nondir_sym_fit$opt$aicc, 
         "nondir_asym" = nondir_asym_fit$opt$aicc,
         "dir_sym" = dir_sym_fit$opt$aicc, 
         "dir_asym" = dir_asym_fit$opt$aicc 
)

### parameter numbers
k = c(
 "er" = 1,
 "dir_sym" = 2, 
 "dir_asym" = 4,
 "nondir_sym" = 3,
 "nondir_asym" = 6
 )

### delta aicc
daicc = sort(aicc  - min(aicc))

### lowest daicc
fir_model = names(daicc[1])
sec_model = names(daicc[2])

### chossing best transition model
if (daicc[2] >= 2) {
  best_model = fir_model
} 
if (daicc[2] < 2) {
  if(k[fir_model] < k[sec_model]){
    best_model = fir_model
  }
  if(k[fir_model] > k[sec_model]){
    best_model = sec_model
  }
}

### infer simmaps
all_maps = phytools::make.simmap(tree = mcc_phylo, 
                                 x = spp_states, 
                                 model = my_models[[best_model]],
                                 pi = c(0.25,0.5,0.25),
                                 nsim = 100
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

tiff("4_graphics/simmap.tiff", units="cm", width=7, height=14, res=600)
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
dev.off()

############################### OUWIE ESTIMATES ################################

### choose a trait!
trait_name = "plant_hei"
### import path
imp_path = paste0("2_trait_results/OUWIE/", trait_name)
### importing model fit
all_best_models = readRDS( paste0(imp_path, "/all_best_models.RDS"))
### importing model parameters
all_best_estimates =  readRDS( paste0(imp_path, "/all_best_estimates.RDS"))

### most frequent best fit 
count_models = table(all_best_models$model)
most_freq_model = names(count_models)[count_models == max(count_models)]
### vector w position of best fit
i_best = which(all_best_models$model == most_freq_model)
### select parameter estimates from most freq best model
freq_best_estimates = all_best_estimates[i_best]
### from list to df
est_df = data.frame( t( sapply(freq_best_estimates,c) ) )

### pick one paramenter
param_name = "theta"
col_names = colnames(est_df)[grepl(param_name, colnames(est_df))]
est = est_df[col_names]

### organzing df
param_df = est %>% 
  pivot_longer(cols = any_of(col_names),
               names_to = "state",
               values_to = "estimate") %>% 
  mutate(
    state = case_when(
      grepl("_1", state) ~ "forest_specialist",
      grepl("_2", state) ~ "generalist",
      grepl("_3", state) ~ "open_specialist"
    )
  )

### outlier boundaries
med = median(param_df$estimate, na.rm = T) 
bound= IQR(param_df$estimate, na.rm = T)*1.5
up_bound = med + bound
dw_bound = med - bound

### cleaning data
clean_param_df = param_df %>% 
  filter(estimate < up_bound & estimate > dw_bound)

### describe
clean_param_df %>% 
  group_by(state) %>% 
  reframe(x = mean(estimate), SD = sd(estimate))

### y axis name
if(trait_name == "height"){
  y_axis_name = "plant height (m)"
}
if(trait_name == "leaf_sla"){
  y_axis_name = "SLA (mm2/mg)"
}
if(trait_name == "seed"){
  y_axis_name = "seed mass (mg)"
}

### parameter symbol
if(param_name == "sigma"){
  param_symbol = "\u03c3" 
}
if(param_name == "theta"){
  param_symbol = "\u03b8" 
}

# plot parameter
plot_param = ggplot(data= clean_param_df, 
                    aes(x=state, 
                        y=estimate, 
                        fill=state)) +
  
  geom_point(aes(color=state),
             position = position_jitter(width = 0.07), 
             size = 0.6, 
             alpha = 0.65) +
  
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25)+

  scale_fill_manual(values=state_cols)+
  scale_colour_manual(values=state_cols)+
  scale_x_discrete(labels= c(
    "rainforest",
    "generalist",
    "open-vegetation"
  )
  ) +
  
  xlab("habitat type") + 
  ylab(paste(param_symbol, y_axis_name, sep = " ") )+
  
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"), 
        axis.title=element_text(size=10,face="bold"), 
        axis.text=element_text(size=8), 
        legend.position = "none") 

### exporting directory
dir_out = "4_graphics/"
file_name = paste(trait_name, param_name, sep ="_")

### export plot
tiff(paste0(dir_out,file_name, ".tiff"), 
     units="cm", width=7, height=6.5, res=600)
print(plot_param)
dev.off()
  
########################### CORRELATION ESTIMATES #############################

### trait name
t1 = "seed_mass"

### directory name
dir_name = paste0("3_trait_results/EVOLVCV/",t1)

### best model list and parameters
fir_model_cor = readRDS(paste0(dir_name,"/fir_model_cor.RDS") )
### best model list and parameters
sec_model_cor = readRDS(paste0(dir_name,"/sec_model_cor.RDS") )

### processing model 1
fir_model_cor = data.frame( cbind("overall", fir_model_cor) )
colnames(fir_model_cor) = c("habitat_range", "cor_values")
fir_model_cor$cor_values = as.numeric(fir_model_cor$cor_values)

### processing model 2
sec_model_cor = sec_model_cor %>% 
  pivot_longer(cols = c(forest_specialist, generalist, open_specialist),
               names_to = "habitat_range",
               values_to = "cor_values"
  )

### plot model 1
model1_plot = ggplot(data = fir_model_cor, 
                     aes(x = cor_values)
  ) +
  
  geom_histogram(position = "identity",
                 fill ="lightgray",
                 color= "black",
                 alpha= 0.5) +
  xlim(c(-1,1)) +
  xlab("Correlation seed mass - SLA") +
  ylab("Frequency")+
  
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"), 
        axis.title=element_text(size=10,face="bold"), 
        axis.text=element_text(size=8), 
        legend.position = c(0.8, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.5,"line")
  ) 

### exporting directory
dir_out = "4_graphics/"
### export plot
tiff(paste0(dir_out,"model_1", ".tiff"), 
     units="cm", width=7.5, height=6.5, res=600)
print(model1_plot)
dev.off()


### plot model 2
model2_plot = ggplot(data = sec_model_cor, 
       aes(x = cor_values, 
           fill = habitat_range)
  ) +
  
  geom_histogram(position = "identity",
                 color= "black",
                 alpha= 0.5) +
  
  scale_fill_manual(values=state_cols,
                    labels= c(
                      "rainforest",
                      "generalist",
                      "open-vegetation"
                    ))+
  
  xlab("Correlation seed mass - SLA") +
  ylab("Frequency")+
  
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"), 
        axis.title=element_text(size=10,face="bold"), 
        axis.text=element_text(size=8), 
        legend.position = c(0.8, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.5,"line")
        ) 

### exporting directory
dir_out = "4_graphics/"
### export plot
tiff(paste0(dir_out,"model_2", ".tiff"), 
     units="cm", width=7.5, height=6.5, res=600)
print(model2_plot)
dev.off()
