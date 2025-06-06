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

### loading trait data
flower_mtx = read.table("0_data/flower_trait_matrix.csv", 
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

################################ trait data ####################################

### flower size per species
flower = flower_mtx %>% 
  group_by(species) %>% 
  reframe(
    flower_size = median(style_height, na.rm=T)
  )
  
### trait values per species
spp_traits = trait_mtx %>% 
  group_by(species) %>% 
  reframe(
    plant_size =  median(plant_height_m, na.rm=T),
    inflor_size = median(inflorescence_length_cm, na.rm=T),
    inflor_rel_size = 100*(inflor_size/(plant_size*100)),
    fruit_size =  median(fruit_height_mm, na.rm=T),
    seed_num = median(seed_number, na.rm=T),
    seed_size = median(seed_height_mm, na.rm=T),
    seed_mass = median(fruit_weight_mg/seed_number, na.rm=T),
    n = n()
  ) %>% 
  left_join(flower, by = "species") %>% 
  left_join(habitat_range, by = "species")

########################## species trait by range #########################

## choose a trait
trait_name = "fruit_size"

### summary
spp_traits %>% 
  group_by(range) %>% 
  reframe(x = mean(fruit_size, na.rm = T),
          stdv = sd(fruit_size, na.rm = T))

summary(aov(fruit_size ~ range, data = spp_traits))

### graphical param
## y axis name
if(trait_name == "plant_size"){
  y_axis_name = "Plant size (m)"
}
if(trait_name == "inflor_size"){
  y_axis_name = "Inflorescence size (cm)"
}
if(trait_name == "flower_size"){
  y_axis_name = "Flower size (mm)"
}
if(trait_name == "fruit_size"){
  y_axis_name = "Fruit size (mm)"
}
if(trait_name == "seed_num"){
  y_axis_name = "Number of seeds per fruit"
}
if(trait_name == "seed_size"){
  y_axis_name = "Seed size (mm)"
}

### plot
trait_plot = ggplot(data= spp_traits,
       aes(x= spp_traits[["range"]],
           y= spp_traits[[trait_name]],
           color = spp_traits[["range"]],
           fill= spp_traits[["range"]])
) +
  
  geom_point(position = position_jitter(width = 0.10), 
             size = 1, 
             alpha = 0.75
  ) +
  
  geom_boxplot(weight = 0.5/length(unique(spp_traits[["range"]])),
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

######################## specimen correlation plots #######################

### traits per sample
samp_traits = trait_mtx %>% 
  reframe(
    species = species,
    plant_size =  plant_height_m,
    inflor_size = inflorescence_length_cm,
    fruit_size =  fruit_height_mm,
    seed_num = seed_number,
    seed_mass = fruit_weight_mg/seed_number,
  ) %>% 
  drop_na() %>% 
  left_join(habitat_range, by = "species")
 
### choose a trait
pred_name = "fruit_size"
resp_name = "seed_mass"

### lm test
summary(lm(log(samp_traits[[resp_name]]) ~ samp_traits[[pred_name]]))

### graphical param
## x axis name
if(pred_name == "plant_size"){
  x_axis_name = "Plant size (m)"
}
if(pred_name == "inflor_size"){
  x_axis_name = "Inflorescence size (cm)"
}
if(pred_name == "flower_size"){
  x_axis_name = "Flower size (mm)"
}
if(pred_name == "fruit_size"){
  x_axis_name = "Fruit size (mm)"
}
if(pred_name == "seed_num"){
  x_axis_name = "Number of seeds per fruit"
}
if(pred_name == "seed_mass"){
  x_axis_name = "Seed mass (mg)"
}

## y axis name
if(resp_name == "plant_size"){
  y_axis_name = "Plant size (m)"
}
if(resp_name == "inflor_size"){
  y_axis_name = "Inflorescence size (cm)"
}
if(resp_name == "flower_size"){
  y_axis_name = "Flower size (mm)"
}
if(resp_name == "fruit_size"){
  y_axis_name = "Fruit size (mm)"
}
if(resp_name == "seed_num"){
  y_axis_name = "Number of seeds per fruit"
}
if(resp_name == "seed_mass"){
  y_axis_name = "Seed mass (mg)"
}

### plot
corr_plot = ggplot(data= spp_traits,
                   aes(x= spp_traits[[pred_name]],
                       y= spp_traits[[resp_name]]) ,
                   )+
                     
  geom_point(aes(color = range),
             size = 2, 
             alpha = 0.75
  ) +
  
  #geom_smooth(method= "lm" , se = F, color = "black")+
  
  xlab(x_axis_name)+ ylab(y_axis_name)+
  
  scale_fill_manual(values=state_cols)+
  scale_colour_manual(values=state_cols)+

  theme(panel.background=element_rect(fill="white"), 
      panel.grid=element_line(colour=NULL),
      panel.border=element_rect(fill=NA,colour="black"),
      axis.title=element_text(size=8,face="bold"),
      axis.text=element_text(size=6),
      legend.position = "none")

## file name
file_name = paste0("4_graphics/corr_",resp_name,"_",pred_name, ".tiff")
### plotting 
tiff(file_name, 
     units="cm", width=7, height= 6.5, res=600)
print(corr_plot)
dev.off()

