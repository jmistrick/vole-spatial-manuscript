# Run using R version 4.3.2 "Eye Holes"

# 04 - Efficacy of deworming treatment
# Author: Janine Mistrick
# Code developed with input from JASMINE S.M. VEITCH & SARAH A. BUDISCHAK
# Associated Publication:
# Effects of Food Supplementation and Helminth Removal on Space Use and Spatial Overlap
# in Wild Rodent Populations.  
  # Mistrick, Janine; Veitch, Jasmine; Kitchen, Shannon; Clague, Samuel; Newman, Brent; 
  # Hall, Richard; Budischak, Sarah; Forbes, Kristian; Craft, Meggan


#Effect of deworm treatment on helminth presence and helminth infection intensity

#Testing on FULL 2021 FEC dataset (1025 entries)
#This code runs with Fecal Egg Count (FEC) data from Dr. Sarah Budischak "prelim_FEC_data.csv"

#load packages
library(here) #v1.0.1
library(tidyverse) #2.0.0
library(lme4) #v1.1-35.1
library(lmerTest) #like lme4 but with pvalues v3.1-3

# #clear environment
rm(list = ls())


############################## LOAD DATA #######################################

#metadata for FEC data
metadata <- readRDS(here("fulltrap21_JAEfinal_withsampID.rds")) %>% 
  group_by(samp_id) %>% slice(1) %>% #remove duplicate entries per occasion (ie week recaps)
  dplyr::select(c("site", "trt", "food_trt", "helm_trt",
                  "season", "month", "occasion",
                  "tag", "samp_id", "sex", "season_breeder"))

#preliminary FEC data from Sarah v.Dec 19 2023
FECdata <- read_csv(here("prelim_FEC_data.csv"))

#pull 2021 FEC data, combine with metadata, keep only complete entries
#helm_trt, pre_post, and nematode.y.n are all factors
all2021FEC <- FECdata %>%
  filter("exclude from FEC analysis" != 1) %>% #remove samples with feces weight less than 0.005
  rename("samp_id" = "Sample ID",
         "nematode.epg" = "nematode epg") %>%
  select(c("samp_id", "nematode.epg", "nematode.y.n")) %>%
  inner_join(metadata, by="samp_id") %>%
  drop_na(nematode.epg) %>% drop_na(nematode.y.n) %>% #remove entries missing FEC data
  mutate(nematode.epg = as.numeric(nematode.epg),
       nematode.y.n = as.factor(nematode.y.n)) %>%
  #remove three duplicate entries (details to follow)
  filter(samp_id != 274) %>% #vole 226170 sampled twice in July, remove second sample ID (eggs detected in both samples)
  filter(samp_id != 340) %>% #vole 219980 sampled twice in Aug, remove first sample ID (no eggs in first sample, eggs in second)
  filter(samp_id != 712) %>% #vole 219809 was sampled twice in Sept, remove first sample ID (eggs detected in both samples)
  filter(samp_id != 64) %>% #remove negative FEC value
  group_by(tag) %>% arrange(occasion, .by_group = TRUE) %>%
  mutate(capt_nbr = row_number()) %>% #add column for the capture number (ie the first, second, third capture of given vole)
  mutate(pre_post = as.factor(case_when(capt_nbr == 1 ~ "pre",
                                        capt_nbr > 1 ~ "post"))) %>%
  # drop_na(sex) %>% drop_na(season_breeder) %>% #model has 10ish voles with sex or repro=NA
  ungroup()

################################################################################################
############################ RUN PREVALENCE & INTENSITY MODELS #################################
################################################################################################
############## here we're running two models, prevalence and intensity #########################
########### models have both pre- and post- measurements, interaction term #####################
################################################################################################

#this code runs on the 'all2021FEC' dataset (models combine pre- and post- data with an interaction by helm_trt) ...
  #using 'helm_trt' and 'pre_post' (both FACTORS) as predictors

#p = prevalence model
#if needed: adjust the optimizer: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
p.prepost <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                   data = all2021FEC, family = binomial) #n=1035
summary(p.prepost) #trt*prepost b=-0.671 (log odds) p=0.0215*
exp(coef(summary(p.prepost))[,1]) #trt*prepost b=0.511 (odds ratio)

# ###### GLMM model diagnostics ######
# #https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
# library(DHARMa)
# #calculate residuals (then run diagnostics on these)
# simulationOutput <- simulateResiduals(fittedModel = p.prepost)
# #qq plot and residual vs fitted
# plot(simulationOutput) #qq look good, scaled residuals look good
# #formal test for overdispersion
# testDispersion(simulationOutput) #looks good, pvalue is large
# #formal test for zero inflation (common type of overdispersion)
# testZeroInflation(simulationOutput) #looks good, pvalue is large
# plotResiduals(simulationOutput, all2021FEC$helm_trt) #looks good
# plotResiduals(simulationOutput, all2021FEC$pre_post) #looks good
# #### OVERALL: DHARMa looks good #####

#------------------------------------------------

#i = intensity (only samples with nematode.epg >0)
i.prepost <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                  data = subset(all2021FEC, nematode.epg > 0)) #n=440
summary(i.prepost) #trt*prepost b=-0.740 p=0.0053**

# #model diagnostics
# plot(i.prepost) #okay... but some patterning
# shapiro.test(resid(i.prepost)) #p=0.00985 pretty small, some issues of non-normality
#
# ######## intensity LMM model diagnostics ############
# #LMM model diagnostics > https://goodekat.github.io/redres/
# # devtools::install_github("goodekat/redres")
# library(redres)
# # creates a plot of the conditional studentized residuals versus the fitted values
# plot_redres(i.prepost, type="std_cond")
# plot_redres(i.prepost, type="pearson_cond") #Pearson conditional residuals
# # creates a residual quantile plot for the error term
# plot_resqq(i.prepost) #looks good
# # creates normal quantile plots for each random effect
# plot_ranef(i.prepost) #also looks good
# ###### OVERALL: model diagnostics look good

#------------------------------------------------

#summary table of mixed model output
#https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
library(sjPlot)
#save it https://stackoverflow.com/questions/67280933/how-to-save-the-output-of-tab-model
tab_model(p.prepost, file="SUPP_Table_S3.doc")
tab_model(i.prepost, file="SUPP_Table S4.doc")


#################################### VISUALIZE MODEL OUTPUT ############################################

## visualize the interaction effects

library(visreg) #https://pbreheny.github.io/visreg/articles/web/overlay.html
library(cowplot)
library(ggtext)

### PREVALENCE MODEL ###

#change pre_post to numeric to get a regression line in visreg output
all2021FEC.num <- all2021FEC %>%
  mutate(pre_post = as.numeric(pre_post)) %>%
  mutate(nematode.y.n_num = case_when(nematode.y.n == "0" ~ 0,
                                      nematode.y.n == "1" ~ 1)) #nematode.y.n from factor to numeric to plot points on continous y axis
p.prepost.num <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                           data = all2021FEC.num, family = binomial)

#visualize infection probability, pretty
png(filename="SUPP_Fig_S3-A.png", height=4, width=6, units="in", res=600)
visreg(p.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
                gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(breaks=c(1,2),
                     expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=11),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="Likelihood of Helminth Infection") +
  annotate(geom = "text", x=1.35, y=1.05, size = 4,
           label = paste("Trt Grp*Trt Stage")) +
  annotate(geom = "text", x=1.35, y=0.97, size = 4,
           label = paste("OR =",
                         round( exp(coef(summary(p.prepost.num))[4,1]), digits=3) )) +
  annotate(geom = "text", x=1.35, y=.90, size = 4,
           label = paste("p =",
                         round( coef(summary(p.prepost.num))[4,4], digits=3) )) +
  geom_jitter(aes(x=pre_post, y=nematode.y.n_num, color=helm_trt),
              data=all2021FEC.num,
              width=0.1, height=0.05,
              size=3, alpha=0.2, shape=16)
dev.off()


### INTENSITY MODEL ###

#pre_post as numeric to look like a regression line
#visreg doesn't like it when the model has a subset command in it #https://github.com/pbreheny/visreg/issues/99
witheggs2021FEC.num <- subset(all2021FEC.num, nematode.epg > 0)
i.prepost.num <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                  data = witheggs2021FEC.num)

#visualize infection intensity, pretty
png(filename="SUPP_Fig_S3-B.png", height=4, width=6, units="in", res=600)
visreg(i.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(breaks=c(1,2),
                     expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="ln(eggs per gram of feces)") +
  annotate(geom = "text", x=1.35, y=8.6, size = 4,
           label = paste("Trt Grp*Trt Stage")) +
  annotate(geom = "text", x=1.35, y=8.0, size = 4,
           label = paste("\u03B2 =",
                         round( coef(summary(i.prepost.num))[4,1], digits=3) )) +
  annotate(geom = "text", x=1.35, y=7.5, size = 4,
           label = paste("p =",
                         round( coef(summary(i.prepost.num))[4,5], digits=3) )) +
  geom_jitter(data=witheggs2021FEC.num, aes(x=pre_post, y=log(nematode.epg), color=helm_trt),
              width=0.05, size=3, alpha=0.2, shape=16)
dev.off()
