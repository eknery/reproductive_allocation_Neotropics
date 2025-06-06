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
  ) %>% 
  left_join(habitat_range, by = "species")

## state colors 
state_cols = c("forestgreen", "brown","darkorange") 
names(state_cols) = levels(habitat_range$range)

################################## PGLS 1 ######################################

### choosing a trait
pred_name = "seed_mass"
pred = spp_traits[[pred_name]]
names(pred) = spp_traits$species

resp_name = "sla"
resp = log(spp_traits[[resp_name]])
names(resp) = spp_traits$species

### fitting models
fit_bm = fitContinuous(phy= mcc_phylo, dat = resp,  model="BM")
fit_ou = fitContinuous(phy= mcc_phylo, dat = resp,  model="OU")

### choosing model aicc
if(fit_bm$opt$aicc < fit_ou$opt$aicc){
  sigsq = fit_bm$opt$sigsq
  spp = names(resp)
  cor_str = corBrownian(sigsq, phy = mcc_phylo, form= ~spp)
}
if(fit_bm$opt$aicc > fit_ou$opt$aicc &
   (fit_bm$opt$aicc - fit_ou$opt$aicc) >= 2 ){
  alpha = fit_ou$opt$alpha
  spp = names(resp)
  cor_str = corMartins(alpha, phy = mcc_phylo, form= ~spp)
}

### fitting pgls
fit_gls1 = gls(resp ~ pred ,
              correlation= cor_str, 
              method = "REML")

summary(fit_gls1)
plot(fit_gls1)

### checking residuals
res1 = resid(fit_gls1)[1:nrow(spp_traits)]
hist(res1)
shapiro.test(res1)

### getting coefficients
intercept1 = fit_gls1$coefficients[["(Intercept)"]]
slope1 = fit_gls1$coefficients[["pred"]]

### plotting
pgls1 = ggplot(data= spp_traits, 
               aes(x= spp_traits[[pred_name]], 
                   y= log(spp_traits[[resp_name]])
                   ) 
       ) +
  
  geom_point(size = 2, alpha = 0.25) + 
  
  geom_abline(intercept = intercept1, 
              slope = slope1, 
              color="black",  
              size= 1,
              linetype="dashed"
              ) +
  
  xlab("seed mass (mg)") + ylab("ln SLA (mm2/mg)")+
  
  theme(panel.background=element_rect(fill="white"), 
      panel.grid=element_line(colour=NULL),
      panel.border=element_rect(fill=NA,colour="black"), 
      axis.title=element_text(size=10,face="bold"), 
      axis.text=element_text(size=8), 
      legend.position = "none") 

### exporting directory
dir_out = "4_graphics/"
tiff(paste0(dir_out,"pgls1", ".tiff"), 
     units="cm", width=7.5, height=6.5, res=600)
print(pgls1)
dev.off()

################################## PGLS 2 ######################################

### fitting pgls
fit_gls2 = gls(resp ~ pred:spp_states ,
               correlation= cor_str, 
               method = "REML"
              )

summary(fit_gls2)
plot(fit_gls2)

### checking residuals
res2 = resid(fit_gls2)[1:nrow(spp_traits)]
hist(res2)
shapiro.test(res2)

### getting coefficients
intercept2 = fit_gls2$coefficients[["(Intercept)"]]
slope2a = fit_gls2$coefficients[["pred:spp_statesforest_specialist"]]
slope2b = fit_gls2$coefficients[["pred:spp_statesgeneralist"]]
slope2c = fit_gls2$coefficients[["pred:spp_statesopen_specialist"]]

### plotting
pgls2 = ggplot(data= spp_traits, 
              aes(x= spp_traits[[pred_name]], 
                  y= log(spp_traits[[resp_name]])
                  ) 
  ) +
  
  geom_point(
    aes(color = range),
    size = 2, 
    alpha = 0.25
    ) + 
  
  scale_color_manual(
    values=state_cols,
    labels= c(
      "rainforest",
      "generalist",
      "open-vegetation"
      )
    )+
  
  geom_abline(
    intercept = intercept2, 
    slope = slope2a, 
    color="forestgreen" ,  
    size= 1,
    linetype="dashed"
  ) +
  
  geom_abline(
    intercept = intercept2, 
    slope = slope2b, 
    color="brown",  
    size= 1,
    linetype="solid"
  ) +
  
  geom_abline(
    intercept = intercept2, 
    slope = slope2c, 
    color="darkorange",  
    size= 1,
    linetype="solid"
  ) +
  
  xlab("seed mass (mg)") + ylab("ln SLA (mm2/mg)")+
  
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"), 
        axis.title=element_text(size=10,face="bold"), 
        axis.text=element_text(size=8), 
        legend.position = "none") 

### exporting directory
dir_out = "4_graphics/"
tiff(paste0(dir_out,"pgls2", ".tiff"), 
     units="cm", width=7.5, height=6.5, res=600)
print(pgls2)
dev.off()

