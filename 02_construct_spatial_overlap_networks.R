# Run using R version 4.3.2 "Eye Holes"

# 02 - Construct Spatial Overlap Networks
# Author: Janine Mistrick
# Associated Publication:
# Effects of Food Supplementation and Helminth Removal on Space Use and Spatial Overlap
# in Wild Rodent Populations.  
  # Mistrick, Janine; Veitch, Jasmine; Kitchen, Shannon; Clague, Samuel; Newman, Brent; 
  # Hall, Richard; Budischak, Sarah; Forbes, Kristian; Craft, Meggan

# This code sources the functions in 02-1_functions_construct_spatial_overlap_networks
# to estimate parameters describing bank vole space use, construct spatial overlap networks representing
# populations of bank voles at a given study site, and calculate network metrics from these networks


# load packages
library(here) #v1.0.1
library(tidyverse) #v2.0.0
library(igraph) #v1.6.0
library(lubridate) #v1.9.3
library(janitor) #v2.2.0

#clear environment
rm(list = ls())


#######----------- GENERATE SPACE USE DISTRIBUTION PARAMETERS ----------###############
#######---------------- CREATE SPATIAL OVERLAP NETWORKS ----------------###############
#######------------------- CALCULATE NETWORK METRICS -------------------###############

#load the 2021 capture data from RDS
fulltrap21 <- readRDS(here("fulltrap21_volecapturedata.rds"))

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

#call the functions
source(here("02-1_functions_construct_overlap_networks.R"))

#generate space use distribution parameters
generate_params(data = fulltrap21,
                params_file = "spaceuse_parameters.rds")
#view output file:
# params21 <- readRDS(here("spaceuse_parameters.rds"))


#create spatial overlap networks from capture data and parameters
create_overlap_networks(data = fulltrap21,
                        centroids_file = "monthly_centroids.rds",
                        params_file = "spaceuse_parameters.rds",
                        networks_file = "overlap_networks.rds")
#view output files:
# centroids21 <- readRDS(here("monthly_centroids.rds"))
# overlapnets21 <- readRDS(here("overlap_networks.rds"))


#calculate network metrics (weighted degree)
calculate_network_metrics(data=fulltrap21,
                          networks_file = "overlap_networks.rds",
                          netmets_file = "network_metrics.rds")
#view output file:
# netmets21 <- readRDS(here("network_metrics.rds"))
