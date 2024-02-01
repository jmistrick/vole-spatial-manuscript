# Run using R version 4.3.2 "Eye Holes"

# 03 - Data analysis, tables, and figures
# Author: Janine Mistrick

#Data analysis and visualization for the main text and supplement of associated publication:

#Effects of Food Supplementation and Helminth Removal on Space Use and Spatial Overlap
# in Wild Rodent Populations.  
  # Mistrick, Janine; Veitch, Jasmine; Kitchen, Shannon; Clague, Samuel; Newman, Brent; 
  # Hall, Richard; Budischak, Sarah; Forbes, Kristian; Craft, Meggan

# load packages
library(here) #v1.0.1
library(tidyverse) #v2.0.0
library(igraph) #v1.6.0
library(ggridges) #0.5.0
library(ggforce) #for geom_circle in ggplot v0.4.1
library(cowplot) #v1.1.2

#clear environment
rm(list = ls())

#load vole capture data from RDS
ft21 <- readRDS(here("fulltrap21_volecapturedata.rds"))

#alternatively, load capture data from csv file and format data columns
fulltrap <- read.csv(here("fulltrap21_volecapturedata.csv")) %>%
  mutate(year = as.numeric(year),
         month = factor(month, levels=c("june", "july", "aug", "sept", "oct")),
         season = factor(season, levels=c("summer", "fall")),
         occ.sess = as.character(occ.sess),
         occasion = as.numeric(occasion),
         session = as.numeric(session),
         site = as.character(site),
         trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm")),
         food_trt = factor(food_trt, levels=c("unfed", "fed")),
         helm_trt = factor(helm_trt, levels=c("control", "deworm")),
         tag = as.character(tag),
         firstcap = factor(firstcap),
         trap = as.character(trap),
         x = as.numeric(x),
         y = as.numeric(y),
         sex = as.factor(sex),
         season_breeder = factor(season_breeder, levels=c("breeder", "nonbreeder")))

##########################################################
################## MAIN TEXT - RESULTS ###################
##########################################################

############ Capture Data Summary #################

#create caps_per_life column - total number of captures of each unique vole
ft21 <- ft21 %>%
  group_by(tag) %>%
  mutate(caps_per_life = length(tag)) %>%
  ungroup()

n_distinct(ft21$tag) #742 individual voles

#summarize mean captures per vole across all 12 sites
onepertag <- ft21 %>% group_by(tag) %>% slice(1) %>% ungroup()
onepertag %>% summarise(mean = mean(caps_per_life),
            sd = sd(caps_per_life),
            min = min(caps_per_life),
            max = max(caps_per_life))

#voles captured at least twice
multicap <- ft21 %>% filter(caps_per_life >= 2) 
n_distinct(multicap$tag) #440 individuals (59.3%) captured at least 2x

#summarize n captures per season
ft21 %>% group_by(season) %>% 
  summarise(n = length(tag)) #905 in summer, 1107 in autumn


##########################################################
################## DATA VISUALIZATIONS ###################
##########################################################


##-------- Visualize Space Use by Sex, Reproductive Status ----------

######## HOME RANGE DISTRIBUTION ###########
#negative sigmoidal curve calculated following Wanelik & Farine 2022:
    # Wanelik, K.M., Farine, D.R. A new method for characterising shared space use networks using animal trapping data.
    # Behav Ecol Sociobiol 76, 127 (2022). https://doi.org/10.1007/s00265-022-03222-5

#where the declining probability (P) of an individual being detected at a distance (d)
#from the centroid of its home range is given by:

P(d) = 1 / (1 + e^(-a-bd))

#where a describes the overall size of the home range
#b describes the steepness of the edge of the home range
#and d is the logarithmic distance from the centroid

# p <- 0.01 #e.g. probability of detection 1%
# a <- params21[14,2]
# b <- params21[14,3]
#
# (log((1/p)-1) + a) / (-b)

##### ^this equation is what I'm using to calculate e.g. the distance from the centroid with 1% probability of finding the vole,
##just to plot some figure showing approx HR size to visualize differences in space use between treatments


####----------- LOAD DATA -----------------

#HR distribution parameters (a and b params)
params21 <- readRDS(here("spaceuse_parameters.rds"))


####----------- GENERATE FIGURE -----------------

#create labellers for figure
trt.labs <- as_labeller(c("unfed_control" = "Unfed-Control",
                          "unfed_deworm" = "Unfed-Deworm",
                          "fed_control" = "Fed-Control",
                          "fed_deworm" = "Fed-Deworm"))

season.labs <- as_labeller(c("summer" = "Summer",
                             "fall" = "Autumn"))

### visualizing the space use to where there is a 1% chance of detecting the vole ###
#### compare space use by sex, reproductive status, treatment, season #####

png(filename = here("Fig_2.png"), height=7, width = 16, units = "in", res=600)
params21 %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b),
         area = paste( round(pi*(rad_0.01^2)*10, digits = 2), "m\u00B2" )) %>%
  separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
  unite(trt, food_trt, helm_trt) %>%
  mutate(x = case_when(sex=="F" ~ 4,
                       sex=="M" ~ 8.5),
         y = case_when(season_breeder=="breeder" ~ 11,
                       season_breeder=="nonbreeder" ~ 7),
         repro = case_when(season_breeder=="breeder" ~ "Reproductive",
                           season_breeder=="nonbreeder" ~ "Non-Reproductive"),
         y_lab = case_when(season_breeder=="breeder" ~ y,
                           season_breeder=="nonbreeder" & season=="summer" ~ y-1.75,
                           season_breeder=="nonbreeder" & season=="fall" ~ y-2.25)) %>%
  mutate(season = factor(season, levels=c("summer", "fall")),
         trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm")),
         repro = factor(repro, levels=c("Reproductive", "Non-Reproductive"))) %>%
  ggplot() +
  geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=sex, linetype=repro), linewidth=0.8) +
  geom_text(aes(x=x, y=y_lab, label=area), hjust=0.5, vjust=0.5, size=5) +
  scale_fill_manual(values=c("#f282a780", "#00d0ff80")) +
  facet_grid(season~trt, labeller=labeller(trt=trt.labs, season=season.labs)) +
  coord_fixed() +
  labs(fill="Sex:", linetype="Reproductive Status:") +
  theme(legend.position="bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size=19),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(order = 1),
         linetype = guide_legend(order = 2)) #adjusts order of legends to show sex first then repro
dev.off()

################# end ######################



############## VISUALIZE SPATIAL OVERLAP NETWORKS ##############
########## PLOT NETWORKS FOR JUNE-OCT, FOR ALL SITES ###########

#load overlap network data
overlap_network_list <- readRDS(here("overlap_networks.rds"))


for(i in 1:length(overlap_network_list)) {
  
  set.seed(2111994)
  
  png(filename = paste("spatial_overlap_", "network-wt.deg0.05_", 
                       names(overlap_network_list)[[i]], "_2021", ".png", sep = ""),
      width=10 , height=3, units="in", res=600)
  
  par(mfrow = c(1,5))
  
  for(j in 1:length(overlap_network_list[[i]])){
    
    data <- overlap_network_list[[i]][[j]]
    
    #plot adj matrix for network
    g <- graph_from_adjacency_matrix(
      data,
      mode = c("undirected"),
      weighted = TRUE,
      diag = FALSE)
    
    #for thresholded edges (visualization purposes only)
    #remove edges with weight less than 0.05 threshold
    g2 <- delete.edges(g, which(E(g)$weight<=0.05))
    
    #keep the layout the same
    #https://www.kateto.net/wp-content/uploads/2016/01/NetSciX_2016_Workshop.pdf
    l <- layout_with_fr(g2)
    
    #for edges of varying thickness
    plot(g2, vertex.size=5, vertex.label=NA,
         layout=l,
         edge.width = ((E(g2)$weight)*3),
         edge.color = "#545454",
         vertex.color = "black",
         main = paste(names(overlap_network_list[[i]])[j]))
    
  }
  
  dev.off()
  
}

########################### end ################################



########### VISUALIZE SPATIAL OVERLAP NETWORK METRICS ###########

#prepare network metrics data with metadata for visualizations

#pull vole sex and tag data
tagsex <- ft21 %>% group_by(tag) %>% slice(1) %>% select(tag, sex, season_breeder)

#join network metrics with metadata on vole sex, and site treatment data
netmets21 <- readRDS(here("network_metrics.rds")) %>% left_join(read.csv(here('site_trts.csv')), by="site") %>%
  unite(trt, food_trt, helm_trt) %>%
  left_join(tagsex, by="tag") %>%
  relocate(c(site, trt, month, n.node, tag, sex, season_breeder), .before=wt.deg)

#labellers for figures
month.labs <- as_labeller(c("june" = "June", "july" = "July", "aug" = "August", "sept" = "September", "oct" = "October"))
breeder.labs <- as_labeller(c("breeder" = "Reproductive", "nonbreeder" = "Non-Reproductive"))

################################ end ####################################


####################### NETWORK SIZE BY MONTH ###########################

#plot network size by month, color by trt
png(filename = 'Fig_4-A.png', width = 12.5, height = 5, units = "in", res=600)
netmets21 %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"  ))) %>%
  group_by(site, month) %>% slice(1) %>%
  ggplot(aes(x=month, y=n.node, color=trt)) +
  geom_jitter(size=2.5, width=0.05) +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A"),
                     name = "Treatment",
                     labels = c("Unfed - Control", "Unfed - Deworm", "Fed - Control", "Fed - Deworm")) +
  scale_x_discrete(labels=c("June", "July", "August", "September", "October")) +
  labs(x="Month", y="Network Size") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=18, margin = margin(0,10,0,0)),
        axis.title.x = element_blank(),
        axis.text = element_text(size=15))
dev.off()

################################### end ##########################################


##################### INDIVIDUAL VOLE SPATIAL OVERLAP #######################

# individual vole spatial overlap (weighted degree) by sex, breeder status for each treatment
png(filename = 'Fig_4-B.png', width = 12, height = 3.5, units = "in", res=600)
netmets21 %>% group_by(trt, month, season_breeder, sex) %>%
  summarise(mean = mean(wt.deg),
            min = min(wt.deg),
            max = max(wt.deg),
            sd = sd(wt.deg)) %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  filter(season_breeder=="breeder") %>%
  ggplot(aes(y=mean, x=trt, color=sex, ymin=min, ymax=max)) +
  geom_point(position=position_dodge(width=0.3), size=3) +
  geom_errorbar(width=0, position = position_dodge(0.3)) +
  scale_color_manual(values=c("#f282a7", "#00d0ff")) +
  facet_grid(season_breeder ~ month,
             labeller = labeller(season_breeder=breeder.labs, month=month.labs)) +
  labs(y="Weighted Degree") +
  coord_flip() + #flip axes so the error bars are horizontal, not vertical
  theme_half_open() + panel_border() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size=15),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(.4, "lines"),
        panel.background = element_rect(fill = "#FAFAfA"))
dev.off()

#####################################  end  ######################################






#####################################################################################
####################### ANALYSIS AND FIGURES FOR SUPPLEMENT #########################
#####################################################################################


####################### sensitivity analysis of unweighted degree ########################

#unweighted (binary) degree (threshold wt > 0.05)
bin05deg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  ggplot( aes(y=trt, x=bin.05.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 20, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Unweighted Degree (Wt Deg > 0.05)", y="Count",
       fill="Treatment") #title="Histogram of Binary Degree Distribution"

#unweighted (binary) degree (threshold wt > 0.01)
bin01deg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  ggplot( aes(y=trt, x=bin.01.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 20, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Unweighted Degree (Wt Deg > 0.01)", y="Count",
       fill="Treatment") #title="Histogram of Binary Degree Distribution"

#unweighted (binary) degree (threshold wt > 0.005)
bin005deg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  ggplot( aes(y=trt, x=bin.005.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 20, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Unweighted Degree (Wt Deg > 0.005)", y="Count",
       fill="Treatment") #title="Histogram of Binary Degree Distribution"

#unweighted (binary) degree (threshold wt > 0.001)
bin001deg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  ggplot( aes(y=trt, x=bin.001.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Unweighted Degree (Wt Deg > 0.001)", y="Count",
       fill="Treatment") #title="Histogram of Binary Degree Distribution"

library(cowplot)
png(filename = 'SUPP_Fig_S1.png', width = 12, height = 12, units = "in", res=600)
plot_grid(bin05deg, bin01deg, bin005deg, bin001deg, ncol = 1)
dev.off()

####################################### end ########################################


#mean and standard deviation of unweighted degree normalized by network size
netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"))) %>%
  mutate(season = case_when(month=="june" ~ "summer",
                            month=="july" ~ "summer",
                            month=="aug" ~ "summer",
                            month=="sept" ~ "fall",
                            month=="oct" ~ "fall")) %>%
  mutate(norm.bin.deg = bin.01.deg/(n.node-1)) %>%
  group_by(season, trt) %>%
  summarise(mean = mean(norm.bin.deg),
            sd = sd(norm.bin.deg))

############################ end ###################################


############# Weighted Degree, Unweighted Degree, Normalized Unweighted Degree by Treatment ############

#ridgeline plot (density) of weighted degree, sites combined by trt
wtdeg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  ggplot( aes(y=trt, x=wt.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Weighted Degree",
       fill="Treatment") #y="Density", title="Histogram of Weighted Degree Distribution"

#binary degree (threshold wt > 0.01)
bindeg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  ggplot( aes(y=trt, x=bin.01.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Unweighted Degree", y="Count",
       fill="Treatment") #title="Histogram of Binary Degree Distribution"

#'normalized' binary degree (# voles you overlap with/# you could overlap with)
normbindeg <- netmets21 %>%
  mutate(month = factor(month, levels=c("june","july", "aug", "sept", "oct"))) %>%
  mutate(trt = factor(trt, levels=c("fed_deworm", "fed_control", "unfed_deworm", "unfed_control"))) %>%
  mutate(norm.bin.deg = bin.01.deg/(n.node-1)) %>%
  ggplot( aes(y=trt, x=norm.bin.deg,  fill=trt)) +
  stat_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x),
                      alpha=0.6, scale=0.9, rel_min_height=0.01) +
  facet_wrap(~month, labeller=month.labs, ncol=5) +
  scale_fill_manual(values = c("#6A3D9A", "#CAB2D6", "#33A02C", "#B2DF8A"),
                    name = "Treatment",
                    labels = c("Fed - Deworm", "Fed - Control", "Unfed - Deworm", "Unfed - Control")) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position="none",
        legend.title = element_text(size=14),
        legend.text=element_text(size=12),
        axis.title.x=element_text(size=18, margin = margin(t = 10, r = 0, b = 10, l = 0)),
        axis.text.x = element_text(size=11),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank()) +
  labs(x="Normalized Unweighted Degree", y="Count",
       fill="Treatment") #title="Histogram of NORMALIZED Binary Degree Distribution"

png(filename = 'SUPP_FIG_S2-B-C-D.png', width = 12, height = 12, units = "in", res=600)
plot_grid(wtdeg, bindeg, normbindeg, ncol = 1)
dev.off()

######################## end #############################



############# Modeling Network Size by Treatment, Month ######################

library(lme4)

#network size per site per month
net_size <- netmets21 %>%
  group_by(month, site) %>%
  slice(1) %>%
  dplyr::select(site,trt,month,n.node) %>%
  separate(trt, c("food", "worms")) %>%
  rename(Site = site,
         network_size = n.node) %>%
  #make month numeric so it's easier to see change through time
  mutate(Month = case_when(month=="june" ~ 0,
                           month=="july" ~ 1,
                           month=="aug" ~ 2,
                           month=="sept" ~ 3,
                           month=="oct" ~ 4)) %>%
  #make food and worms numeric so the reference group is what we want
  mutate(Food_Addition = case_when(food=="fed" ~ 1,
                                   food=="unfed" ~ 0),
         Helminth_Removal = case_when(worms=="control" ~ 0,
                                      worms=="deworm" ~ 1))

# #run a linear model
# netsize_mod <- lm(n.node ~ food + worms + month , data=net_size)
# summary(netsize_mod)
# plot(netsize_mod) #ehh things might be a little funky
# #formal test!
# shapiro.test(residuals(netsize_mod)) #yeah, residuals aren't normal
#
# #just to check: look at the distribution of the network sizes
# hist(net_size$n.node) #kind of skewed
#
# #what if we log-transform n.node?
# #run a linear model
# netsize_mod <- lm(log(n.node) ~ food + worms + month , data=net_size)
# summary(netsize_mod) #looks better (estimates are scaled 0-1)
# plot(netsize_mod) #looks better
# # but confirm with a formal test!
# shapiro.test(residuals(netsize_mod)) #yaas the pvalue is large, we gucci
#
# #and confirm
# hist(log(net_size$n.node)) #yep, more normal shaped

#next step, see if we can make it run in a lmer
netsize_mod <- lmer(log(network_size) ~ Food_Addition + Helminth_Removal + Month + (1|Site), data=net_size)
summary(netsize_mod)

#summary table of mixed model output
#https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
library(sjPlot)
tab_model(netsize_mod, file="SUPP_Table_S2.doc")

# ## MODEL DIAGNOSTICS - CODE FROM Jasmine Veitch (this does the same as below but baseR)
# plot(netsize_mod)
# plot(resid(netsize_mod))
# qqnorm(resid(netsize_mod))
# qqline(resid(netsize_mod))
# hist(resid(netsize_mod))
# #construct CI
# confint(netsize_mod)
# #chi square for pvalues (Ben Bolker said to do this)
# drop1(netsize_mod, test="Chisq")
# #easy code for figures
# library(visreg)
# visreg(netsize_mod) #plotting partial residuals, change one parameter, control the others
# #shows the figures for each of the predictors
# #shows the predicted values based on the model
# ##use this after running the model (look at the model visually in addition to the model summary)
# 
# 
# ## MODEL DIAGNOSTICS - https://dfzljdn9uc3pi.cloudfront.net/2020/9522/1/MixedModelDiagnostics.html <- #MattSilkfortheWin
# #STEP 1: plot standardized residuals vs fitted values
# #resid vs fitted check for homoscedasticity
# plot(netsize_mod, resid(., scaled=TRUE) ~ fitted(.),
#      abline = 0,pch=16,xlab="Fitted values",ylab="Standardised residuals") #looks like a nice cloud
# 
# #STEP 1.5: plot residuals vs explanatory values, tests for linearity
# #In our case we plot standardised residuals against food
# plot(netsize_mod, resid(., scaled=TRUE) ~ Food_Addition,
#      abline = 0, pch=16, xlab="Food",ylab="Standardised residuals")
# #and plot standardised residuals against worms
# plot(netsize_mod, resid(., scaled=TRUE) ~ Helminth_Removal,
#      abline = 0, pch=16, xlab="Worms",ylab="Standardised residuals")
# 
# #STEP 2: Residuals split by random effect groupings (check for homosce... in RE group too)
# #We can use a default plot of the standardised residuals versus fitted
# # #But add | <random effect> to see it per level of the random effect
# # plot(netsize_mod, resid(., scaled=TRUE) ~ fitted(.)| site,
# #      abline = 0,pch=16,xlab="Fitted values",ylab="Standardised residuals")
# #   #hard to say since only 5 points per group
# #You can also compare the difference in residuals between random effect levels using default boxplots
# plot(netsize_mod, as.factor(Site) ~ resid(., scaled=TRUE),
#      abline=0,pch=16,xlab="Standardised residuals",ylab="Sites")
# 
# #STEP 3: QQ plot to make sure residuals are normally distributed
# #qq plot - check that residuals are normally distributed
# qqnorm(resid(netsize_mod),pch=16)
# qqline(resid(netsize_mod)) #meh, decent?
# 
# #STEP 4: Leverage and Cook's Distance - check for INFLUENTIAL DATA POINTS
# #Calculate leverage
# lev<-hat(model.matrix(netsize_mod))
# #Plot leverage against standardised residuals
# plot(resid(netsize_mod,type="pearson")~lev,las=1,ylab="Standardised residuals",xlab="Leverage")
# #Calculate Cook's Distance
# cd<-cooks.distance(netsize_mod)
# max(cd) #If Cook's distance is greater than 1 this highlights problematic datapoints
# #Plot leverage and Cook's distance together
# par(mfrow=c(1,1))
# plot(lev,pch=16,col="red",ylim=c(0,1.2),las=1,ylab="Leverage/Cook's distance value")
# points(cd,pch=17,col="blue")
# points(x=150,y=1.1,pch=16,col="red")
# points(x=150,y=0.9,pch=17,col="blue")
# text(x=155,y=1.1,"Leverage",adj=c(0,0.5))
# text(x=155,y=0.9,"Cook's distance",adj=c(0,0.5))
# 
# #STEP 5: Inspect the RE distribution (should be Gaussian distributed)
# #It is possible to do this using a histogram but this can be unclear if there are few random effect levels as in our example
# #A dotplot is a more effective way of doing this.
# lattice::dotplot(ranef(netsize_mod, condVar=TRUE))

##################################### end #########################################
