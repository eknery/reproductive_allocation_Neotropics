### require packages
library(tidyverse)
library(PupillometryR)
library(ggpubr)
library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)

# my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")

############################ plotting extent species data #########################

### loading data
spp_rao = read.table("2_environmental_heterogeneity/spp_rao.csv", header =T, sep=",",  na.strings = "NA", fill=T)
spp_altitude = read.table("2_environmental_heterogeneity/spp_altitude.csv", header =T, sep=",",  na.strings = "NA", fill=T)
spp_hvolumes = read.table("3_hypervolume_inference/spp_hvolumes.csv", header =T, sep=",",  na.strings = "NA", fill=T)
spp_range = read.table("4_range_inference/spp_range.csv", header =T, sep=",",  na.strings = "NA", fill=T)
  
# organizing dataset to plot
spp_dataset = data.frame(spp_rao, spp_altitude$altitude, spp_hvolumes$hvolume, spp_range$range)
colnames(spp_dataset)[4:6] = c("altitude", "hvolume","range")

tiff("6_graphs/qvalues_geographic_distribution.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= spp_dataset, aes(x=state, y=rao, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("median Q value")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")
dev.off()


### overall model
# intercept = -2.3984521 slope= 0.536290

xs = range(spp_dataset$altitude*10^-3)
ys = (range(spp_dataset$altitude*10^-3)*0.536290) + -2.3984521
xys = data.frame(x1 = xs[1], x2 = xs[2], y1 = ys[1], y2 = ys[2])

tiff("6_graphs/log_qvalues_altitude.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= spp_dataset, aes(x=altitude*10^-3, y=log(rao) )) +
  geom_point(aes(color=state),size = 2, alpha = 0.35) +
  # geom_smooth(method="lm", formula = y ~ x, se=T, fill="gray", color="black", alpha=0.50)+
  geom_segment(data = xys, color="black", size=1.2,  aes(x = x1, y = y1, xend = x2, yend = y2))+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("median altitude (km)")+ ylab("log median Q value")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=8),legend.position = "none")
dev.off()

tiff("6_graphs/hvolumes_geographic_distribution.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= spp_dataset, aes(x=state, y=hvolume, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("hypervolume size")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")
dev.off()

### overall model 
# intercept = 2.7522638 slope = 0.4996889
xs =  range(spp_dataset$hvolume)
ys = (range(spp_dataset$hvolume)*0.4996889) + 2.7522638
xys = data.frame(x1 = xs[1], x2 = xs[2], y1 = ys[1], y2 = ys[2])

tiff("6_graphs/cubic_convex_area_hvolumes.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= spp_dataset, aes(x=hvolume, y=range^(1/3) ) ) +
  geom_point(aes(color=state),size = 2, alpha = 0.35) +
  #geom_smooth(method="lm", formula = y ~ x, se=T, aes(fill=state, color=state), alpha=0.35)+
  geom_segment(data = xys, color="black", size=1.2,  aes(x = x1, y = y1, xend = x2, yend = y2))+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("hypervolume size")+ ylab("convex hull area (km2)")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=8),legend.position = "none")
dev.off()

### model per each state
#aF,  intercept = 0.672154, slope= 0.733839
#Afother,  intercept = 3.485769+0.672154, slope= 0.440802
#other, intercept = 3.501537+0.672154, slope= 0.312107

hv_ranges = aggregate(spp_dataset$hvolume, by=list(spp_dataset$state), range)

af_x= hv_ranges$x[1,]
af_y= (hv_ranges$x[1,]*0.733839) + 0.672154
af_xy = data.frame(x1 = af_x[1], x2 = af_x[2], y1 = af_y[1], y2 = af_y[2])

afot_x= hv_ranges$x[2,]
afot_y= (hv_ranges$x[2,]*0.440802) + 3.485769+0.672154
afot_xy = data.frame(x1 = afot_x[1], x2 = afot_x[2], y1 = afot_y[1], y2 = afot_y[2])

ot_x= hv_ranges$x[3,]
ot_y= (hv_ranges$x[3,]*0.312107) + 3.501537+0.672154
ot_xy = data.frame(x1 = ot_x[1], x2 = ot_x[2], y1 = ot_y[1], y2 = ot_y[2])

tiff("6_graphs/cubic_convex_area_hvolumes_state.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= spp_dataset, aes(x=hvolume, y=range^(1/3) ) ) +
  geom_point(aes(color=state),size = 2, alpha = 0.35) +
  #geom_smooth(method="lm", formula = y ~ x, se=T, aes(fill=state, color=state), alpha=0.35)+
  geom_segment(data = af_xy, color="#1E88E5", size=1.2,  aes(x = x1, y = y1, xend = x2, yend = y2))+
  geom_segment(data = afot_xy, color="#FFC107", size=1.2,  aes(x = x1, y = y1, xend = x2, yend = y2))+
  geom_segment(data = ot_xy, color="#D81B60", size=1.2,  aes(x = x1, y = y1, xend = x2, yend = y2))+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("hypervolume size")+ ylab("convex hull area (km2)")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=8),legend.position = "none")
dev.off()


tiff("6_graphs/convex_area_geographic_distribution.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= spp_dataset, aes(x=state, y=sqrt(range), fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic state")+ ylab("sqrt( range size, km2 )")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")
dev.off()

############################### plotting phylogenetic data ############################

library(ape)
library(phytools)

# mcc tree
mcc = read.tree("0_data/mcc_phylo.nwk")
n_tips = Ntip(mcc)
n_inner_nodes = mcc$Nnode

# DEC probabilities
node_probs_mcc = read.table("5_comparative_analysis/DEC_node_probs_mcc.csv",sep=",", h=T)

# DEC states
node_states_mcc = read.table("5_comparative_analysis/DEC_node_states_mcc.csv", sep=",", h=T)

# hypervolume data
spp_hvolumes = read.table("3_hypervolume_inference/spp_hvolumes.csv", header =T, sep=",",  na.strings = "NA", fill=T)
hvolumes = spp_hvolumes$hvolume
names(hvolumes) = spp_hvolumes$species

###range data
spp_range = read.table("4_range_inference/spp_range.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# normal range
ranges = spp_range$range
names(ranges) = spp_range$species

# cubic range
cubic_ranges = spp_range$range^(1/3)
names(cubic_ranges) = spp_range$species

### visual parameters
# tip states
extent_states = node_states_mcc$x[1:n_tips]
names(extent_states) = mcc$tip.label

# tip states probs
tip_states_probs = node_probs_mcc[1:n_tips, ]

# ancestral state probs
inner_node_probs = node_probs_mcc[(1+n_tips):(n_tips+n_inner_nodes),]

# state colors
state_cols=c( "#1E88E5","#D81B60", "#FFC107")
names(state_cols)=c("AF",  "other", "AFother")

# bar colors
bar_cols=extent_states
bar_cols[extent_states == "AF"] = "#1E88E5"
bar_cols[extent_states == "AFother"] = "#FFC107"
bar_cols[extent_states == "other"] = "#D81B60"
names(bar_cols)=mcc$tip.label

#### plotting map
tiff("6_graphs/dec_mcc_ranges.tiff", units="in", width=4, height=6, res=600)
  plotTree.wBars(ladderize(mcc),ranges,type="phylogram",fsize=0.5, col=bar_cols, lmethod="plotTree")
  tiplabels(pie=tip_states_probs, piecol=state_cols, cex=0.3)
  nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes), pie= inner_node_probs, piecol=state_cols, cex=0.8)
  axisPhylo(pos=c(0.5), font=3, cex.axis=0.5)
dev.off()


########################## plotting observed traits & model estimates ###########################

### best-estimates
best_estimates_table = read.table("5_comparative_analysis/hv_best_estimates.csv", sep=",", h=T)

# organizing best-fit estimates
af = best_estimates_table[-46,c(1,2,7,8)]
afot = best_estimates_table[-46,c(3,4,9,10)]
ot = best_estimates_table[-46,c(5,6,11,12)]
colnames(af) = colnames(afot) = colnames(ot) = c("theta", "se", "alpha", "sigma_sq")
state = c( rep("AF", nrow(af)), rep("AFother", nrow(afot)), rep("other", nrow(ot) ) )
model_estimates = data.frame(state, rbind(af, afot, ot))

aggregate(model_estimates[,2:5], by=list(model_estimates[,1]), mean)
aggregate(model_estimates[,2:5], by=list(model_estimates[,1]), sd)

### BM simulations
bm_simulations = read.table("5_comparative_analysis/hv_bm_simulations.csv", sep=",", h=T)

# organizing bm simulations
model = rep("bm1",nrow(bm_simulations))
bm_simulations = data.frame(model, bm_simulations)
colnames(bm_simulations)[2] = "trait"

## observed hvolumes
regimes = read.table("4_range_inference/regimes_ranges.csv", sep=",", h=T)

# summarizing hvolumes by state
means=aggregate(regimes[,3], by=list(regimes[,2]),mean)
sds=aggregate(regimes[,3], by=list(regimes[,2]),sd)
ns=aggregate(regimes[,3], by=list(regimes[,2]),length)
ses=sds[,2]/sqrt(ns[,2])
summary_regimes=data.frame(means,sds[,2],ses)
colnames(summary_regimes)=c("state","mean","sd","se")

theta_plot= ggplot(data= model_estimates, aes(x=state, y=theta, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.20), size = 1, alpha = 0.1) +
  geom_boxplot(width = 0.50, outlier.shape = NA, alpha = 0.25)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("")+ ylab("theta")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other", "other" = "Outside AF"))+
  #ylim(c(0,750))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=6),legend.position = "none")

sigma_plot= ggplot(data= model_estimates, aes(x=state, y=sqrt(sigma_sq), fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.20), size = 1, alpha = 0.1) +
  geom_boxplot(width = 0.50, outlier.shape = NA, alpha = 0.25)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic\ndistribution")+ ylab("sigma")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other", "other" = "Outside AF"))+
  #ylim(c(0,3000))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=6),legend.position = "none")

bm_plot= ggplot(data=bm_simulations,aes(x=trait)) +
  geom_density(alpha=0.5,color="lightgray", fill="lightgray")+
  scale_colour_manual(values=mycols)+
  scale_fill_manual(values=mycols)+
  geom_vline(data=summary_regimes,aes(xintercept = mean),linetype="solid",colour=c(mycols),size=0.75)+
  geom_vline(data=summary_regimes,aes(xintercept = mean-se),linetype="dotted",colour=c(mycols),size=0.45)+
  geom_vline(data=summary_regimes,aes(xintercept = mean+se),linetype="dotted",colour=c(mycols),size=0.45)+
  labs(x="mean \n convex hull area", y="density")+
  #xlim(c(0,25))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=8),legend.position = "none")

tiff("6_graphs/evo_model_estimates_hv.tiff", units="in", width=2.5, height=6, res=600)
  ggarrange(theta_plot,sigma_plot, bm_plot, nrow=3,ncol=1)
dev.off()

############################# plotting ancestral species data ###########################

### my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")

### ancestral data
anc_data = read.table("5_comparative_analysis/anc_range_data_OUMV.csv", sep=",", h=T)

# time boundaries
old_age = min(anc_data$anc_node_ages)
new_age = round(max(anc_data$anc_node_ages), 3)

# dividing time in intervals
intervals = anc_data$anc_node_ages
breaks = seq(new_age, old_age, by= (old_age - new_age)/10)
for (i in 1:length(breaks)){
  intervals[which(anc_data$anc_node_ages > breaks[i+1] & anc_data$anc_node_ages < breaks[i])] = (breaks[i] + breaks[i+1])/2
}
anc_data = data.frame(anc_data, intervals)

# removing outliers
anc_data = anc_data[anc_data$anc_range > 0 & anc_data$anc_range < 270,]

# summarizing range by state across intervals
list_anc_data = split(anc_data, f= anc_data$state)
all_summary_ranges = data.frame()
for (i in 1:length(list_anc_data)){
  central_ranges = aggregate(list_anc_data[[i]]$anc_range, by= list(list_anc_data[[i]]$intervals), mean)
  dispersion_ranges  = aggregate(list_anc_data[[i]]$anc_range, by= list(list_anc_data[[i]]$intervals), function(x){sd(x)/sqrt(length(x)) })
  state = rep(names(list_anc_data)[i], nrow(central_ranges) )
  summary_ranges = cbind(state, central_ranges, dispersion_ranges[,-1])
  all_summary_ranges = rbind(all_summary_ranges, summary_ranges)
}
colnames(all_summary_ranges) = c("state", "age", "central", "dispersion")

# export
write.table(all_summary_ranges, "6_graphs/all_summary_ranges.csv", sep=",", quote=F, row.names=F)

tiff("6_graphs/convex_hull_area_by_time_OUMV.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= anc_data, aes(x=anc_node_ages, y=anc_range, fill=state) ) +
  #geom_point(aes(color=state),size = 1, alpha = 0.10) +
  geom_smooth(se= T, aes(color=state, fill=state), alpha=0.25, size=1)+
  geom_vline(xintercept = -2.58, linetype="dashed",color = "black", size=0.75) +
  scale_colour_manual(values=mycols)+
  scale_fill_manual(values=mycols)+
  xlab("time before present (m.y.a.)")+ ylab("convex hull area")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10), legend.position = "none")
dev.off()

tiff("6_graphs/mean_convex_hull_area_by_time_OUMV.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= all_summary_ranges, aes(x=age, y=central, group= state, color=state) ) +
  geom_point(size = 1, alpha = 1) +
  geom_line(size=1)+
  geom_errorbar(size=0.75, width=0, aes(ymin=central-dispersion, ymax=central+dispersion))+
  geom_vline(xintercept = -2.58, linetype="dashed",color = "black", size=0.75) +
  scale_colour_manual(values=mycols)+
  scale_x_continuous(breaks=seq(-11,0, by=1), labels=as.character(seq(-11,0, by=1)) )+
  ylim(c(50,200))+
  xlab("time before present (m.y.a.)")+ ylab("convex hull area (km2)")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=9), legend.position = "none")
dev.off()


write.table(all_summary_ranges, "6_graphs/all_summary_ranges.csv", sep=',', row.names = F, quote= F)
