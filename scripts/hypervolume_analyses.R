### packages
if (!require("raster")) install.packages("raster"); library("raster")
if (!require("rgeos")) install.packages("rgeos"); library("rgeos")
if (!require("hypervolume")) install.packages("hypervolume"); library("hypervolume")
if (!require("cli")) install.packages("cli"); library("cli")
if (!require("dismo")) install.packages("dismo"); library("dismo")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("geiger")) install.packages("geiger"); library("geiger")
if (!require("nlme")) install.packages("nlme"); library("nlme")

#loading spp coordinates
spp_points=read.table("0_data/spp_points_7km.csv", 
                      header =T, sep=",",  na.strings = "NA", fill=T)

# loading phylogenetic tree
mcc_phylo =read.tree("0_data/mcc_phylo.nwk")
phylo_trees = read.tree("0_data/100_rand_phylos.nwk")

# loading raster layers
ras1 = raster("0_data/rasters/temperature_diurnal_range")
ras2 = raster("0_data/rasters/precipitation_seasonality")
ras3 = raster("0_data/rasters/solar_radiation")
ras4 = raster("0_data/rasters/soil_pH.gri")
env_ras= stack(ras1,ras2,ras3,ras4)

# all species names
all_spp_names = sort(unique(spp_points$species))

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

# define states threshold
high_ths = 0.90
low_ths = (1 - high_ths)
eco_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
eco_states[af_percentage >= high_ths] = "specialist"
eco_states[af_percentage < high_ths] = "generalist"

names(eco_states) = spp_count_domain$species

# sourcing other functions
source("scripts/function_sister_pairs.R")

################################## data preparation ###########################

### croping rasters
#set extent
ext=extent(-100, -33, -34, 24)
allras_crop=crop(env_ras[[1]], ext)
for (i in 2:length(env_ras@layers) ){
  ras_crop =  crop(env_ras[[i]], ext)
  allras_crop= raster::stack(allras_crop, ras_crop)
}
crop_env_ras = allras_crop


### keeping mean and sd values
ras_means= cellStats(crop_env_ras, stat='mean', na.rm=TRUE)
ras_sds = cellStats(crop_env_ras, stat='sd', na.rm=TRUE)

### converting rasters to z-scores
scale_env_ras = crop_env_ras[[c(1:nlayers(env_ras))]]
for (i in 1:nlayers(scale_env_ras)) { 
  scale_env_ras[[i]] = (scale_env_ras[[i]] - cellStats(scale_env_ras[[i]], 'mean')) / cellStats(scale_env_ras[[i]], 'sd') 
}

### treating spp occurrences
# extracting spp z-values
scale_spp_env = raster::extract(scale_env_ras,spp_points[,2:3])
scale_spp_env = data.frame(spp_points$species, scale_spp_env)
colnames(scale_spp_env)[1] = "species"

# removing spp NAs
spp_nas = c()
for (i in 2:length(scale_spp_env[1,])){
  spp_nas= c(spp_nas, which(is.na(scale_spp_env[,i])) )
  }
if (length(spp_nas) > 0){
  scale_spp_env = scale_spp_env[-spp_nas,]
  }

############################ fitting hypervolumes ######################

# loading threshold values
mean_hv_performance=read.table("1_hypervolume_inference/mean_hv_performance.csv", 
                               header =T, sep=",",  na.strings = "NA", fill=T)

# getting best thresholds for each species
best_ths = rep(NA, length(all_spp_names))
for(sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  mean_sp_performance = mean_hv_performance[mean_hv_performance$species==sp_name,]
  sp_th = mean_sp_performance$threshold[mean_sp_performance$tss == max(mean_sp_performance$tss)] 
  best_ths[index] = sp_th[1]
}
names(best_ths) = all_spp_names

# fitting models per species
all_spp_hv = list()
for (sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  sp_data = scale_spp_env[scale_spp_env$species== sp_name, 2:5]
  sp_ths = best_ths[names(best_ths) == sp_name]
  sp_band = hypervolume::estimate_bandwidth(sp_data)
  sp_hv = hypervolume::hypervolume_gaussian(sp_data,  
                                            samples.per.point = ceiling((10^(3 + sqrt(ncol(sp_data))))/nrow(sp_data)), 
                                            kde.bandwidth = sp_band,
                                            sd.count = 3,
                                            chunk.size = 100,
                                            quantile.requested = sp_ths, 
                                            quantile.requested.type = "probability")
  all_spp_hv[[index]] = sp_hv
}
names(all_spp_hv) = all_spp_names

### exporting hypervolumes
saveRDS(all_spp_hv, "1_hypervolume_inference/all_spp_hv.RDS")

############################# sister hv comparison ############################

  # phylogenetic distance
  phylo_distance = cophenetic(mcc_phylo)
  # sister taxa
  sister_taxa_list = sister_pairs(phylo_distance)
  ### sister divergence time
  sister_divergence = c()
  for(focal_sp in all_spp_names){
    focal_sp_dists= round(phylo_distance[focal_sp,],5)
    min_phylo_dist = round( min(phylo_distance[focal_sp,][phylo_distance[focal_sp,] != 0]), 5)
    sister_divergence = c(sister_divergence, min_phylo_dist)
  }
  ### sister hv comparison
  sister_hv_comparison = matrix(0, nrow=length(all_spp_names), ncol=4)
  for (sp_name in all_spp_names){
    index = which(sp_name == all_spp_names)
    # sp hv
    sp_hv = all_spp_hv[[sp_name]]
    # sister hv
    sister_name = sister_taxa_list[[sp_name]]
    if (length(sister_name) > 1){
      one_sister_name = sister_name[1]
      sister_hv = all_spp_hv[[one_sister_name]]
      for (i in 2:length(sister_name) ){
        one_sister_name = sister_name[i]
        one_sister_hv = all_spp_hv[[one_sister_name]]
        sister_set = hypervolume::hypervolume_set(sister_hv, one_sister_hv, check.memory = F)
        sister_hv = sister_set[[4]]
      }
    } else{
      sister_hv = all_spp_hv[[sister_name]]
    }
    # getting sp and sister volumes
    sister_hv_comparison[index,1] = sp_hv@Volume
    sister_hv_comparison[index,2] = sister_hv@Volume
    # getting intersection and union volumes 
    hv_set = hypervolume::hypervolume_set(sp_hv, sister_hv, check.memory = F)
    sister_hv_comparison[index,3] = hv_set[[3]]@Volume # taking intersection
    sister_hv_comparison[index,4] = hv_set[[4]]@Volume # taking union
  }
  # minimal hypervolume per sister pair
  min_hv = apply(sister_hv_comparison[,1:2], MARGIN = 1, FUN= min)
  # dataframe
  sister_hv_comparison = data.frame(all_spp_names, 
                                    sister_hv_comparison, 
                                    min_hv, 
                                    sister_divergence)
  colnames(sister_hv_comparison) =c("species", 
                                    "sp_hv", 
                                    "sister_hv", 
                                    "intersection", 
                                    "union", 
                                    "minimal_hv",
                                    "divergence_time")
  #exporting
  write.table(sister_hv_comparison, 
              "1_hypervolume_inference/sister_hv_comparisons.csv", 
              sep=',', quote=F, row.names=F)
 
  

######################## phylo GLS #########################

### improting sister hv metrics
sister_hv_comparison = read.table("1_hypervolume_inference/sister_hv_comparisons.csv", 
                                    sep=',', h=T)

### calculating niche dissimilarity  
hv_disim = 1 - (sister_hv_comparison$intersection/sister_hv_comparison$union)
names(hv_disim) = sister_hv_comparison$species

### divergence time
div_time = sister_hv_comparison$divergence_time
names(div_time) = sister_hv_comparison$species

### comparing models
fit_wn = fitContinuous(phy= mcc_phylo, hv_disim,  model="white")
fit_bm = fitContinuous(phy= mcc_phylo, hv_disim,  model="BM")
fit_ou = fitContinuous(phy= mcc_phylo, hv_disim,  model="OU")

### aicc
aicc = c(fit_wn$opt$aicc, fit_bm$opt$aicc, fit_ou$opt$aicc) 
names(aicc) = c("WN", "BM", "OU")
aicc

### correlation structure
#cor_obj = corBrownian(phy= mcc_phylo)
cor_obj = corMartins(value = fit_ou$opt$alpha, phy = mcc_phylo, fixed = F)

### null model
fit_gls0 = gls(hv_disim ~ div_time , 
              correlation= cor_obj, 
              method = "ML")

summary(fit_gls0)
plot(fit_gls0)

### alternative model
fit_gls1 = gls(hv_disim ~ div_time / eco_states,  
               correlation= cor_obj, 
               method = "ML")

summary(fit_gls1)
plot(fit_gls1)

### alternative model
fit_gls2 = gls(hv_disim ~ eco_states,  
               correlation= cor_obj, 
               method = "ML")

summary(fit_gls2)
plot(fit_gls2)

## fast plot
boxplot( hv_disim ~ eco_states)

