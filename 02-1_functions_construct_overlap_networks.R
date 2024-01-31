# 02-1 - Functions to construct spatial overlap networks
# Author: Janine Mistrick
# Associated Publication:
# Effects of Food Supplementation and Helminth Removal on Space Use and Spatial Overlap
# in Wild Rodent Populations.  
  # Mistrick, Janine; Veitch, Jasmine; Kitchen, Shannon; Clague, Samuel; Newman, Brent; 
  # Hall, Richard; Budischak, Sarah; Forbes, Kristian; Craft, Meggan

# This code provides the functions necessary to run the code in 02_construct_spatial_overlap_networks
# and to estimate parameters describing bank vole space use, construct spatial overlap networks representing
# populations of bank voles at a given study site, and calculate network metrics from these networks

#File includes the functions:
  # generate_params()
  # create_overlap_networks()
  # calculate_network_metrics()

##-----------------------------------------------------------------

## Function for generating the a and b parameters by season/trt/sex/season_breeder to define the distributions for vole HRs
## Inputs: data = FULL fulltrap dataframe that contains all capture data for a given year
##          params_file = file name and extension in " " for the params_file to be generated
## Output: dataframe of parameters (3 columns for seasontrtsexbreeder (STSB), "a" param, "b" param)

### this function is only saving the a and b parameters - NOT the locations of the seasonal centroids
    ###(because I'm using monthly centroids to construct the networks and save those locations in the construct_overlap_networks function)

generate_params <- function(data, params_file) {

  ##---------------- LOAD THE CAPTURE DATA ----------------------

  #load the fulltrap dataset
  fulltrap <- data

  #pull all the traps and their x, y coordinates, save as df with (trapID, x, y)
  traps <- fulltrap %>% group_by(trap) %>% slice(1) %>%
    select(trap, x, y) %>%
    arrange(trap)

  #save a vector of the site names (alphabetical order)
  site_names <- unique(fulltrap$site) %>% sort()

  #save a vector of trts
  trt <- unique(fulltrap$trt) %>% sort()

  #save a vector of seasons
  season <- unique(fulltrap$season) %>% sort()

  #save a vector of season_breeder
  season_breeder <- unique(fulltrap$season_breeder) %>% sort()

  ##---------------- CREATE A LIST of CAPTURE DATA (nested by season, trt) --------------------

  #a slick little something from stackoverflow to construct a nested list in one go #blessed
  season_breeder_list <- lapply(split(fulltrap, fulltrap$season, drop = TRUE),
                                function(x) split(x, x[["season_breeder"]], drop = TRUE))

  #create new list to hold nested site, month capture
  sbt_list <- list()

  for(i in 1:length(season_breeder_list)){
    # for each season (summer, fall)
    season_list <- list()

    for(j in 1:length(season_breeder_list[[i]])){
      # for breeders and nonbreeders

      trt_list <- split(season_breeder_list[[i]][[j]], season_breeder_list[[i]][[j]]$trt)
      season_list[[j]] <- trt_list
    }

    names(season_list) <- season_breeder
    sbt_list[[i]] <- season_list

  }

  names(sbt_list) <- season

  ### SBT_LIST (season, breeder, treatment=SBT) is a nested list of length 2
  ### 1e level is season (summer/fall) (2)
  ### 2e level is breeding status in that season (breeder / nonbreeder)
  ### 3e level is all the treatments (4)
  ### df contain data for voles from all grids per trt (3) and months per season combined together!!

  ##------------------------------------

# Following code from Wanelik & Farine 2022
  # Wanelik, K.M., Farine, D.R. A new method for characterising shared space use networks using animal trapping data.
  # Behav Ecol Sociobiol 76, 127 (2022). https://doi.org/10.1007/s00265-022-03222-5

params_list <- list()

  for(i in 1:length(sbt_list)){

    print(i)

    params_list[[i]] <- list()

    for(k in 1:length(sbt_list[[i]])){

      print(k)

      params_list[[i]][[k]] <- list()

      for(j in 1:4){
        #1:4 for the four treatments

        print(j)

        data <- sbt_list[[i]][[k]][[j]] %>%
          select(tag, sex, trap, x, y) %>%
          group_by(tag, trap) %>%
          mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
          ungroup()

        season_lab <- unique(sbt_list[[i]][[k]][[j]]$season)
        breeder_lab <- unique(sbt_list[[i]][[k]][[j]]$season_breeder)
        trt_lab <- unique(sbt_list[[i]][[k]][[j]]$trt)

        tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
        tagsex <- data %>% group_by(tag) %>% slice(1) %>% select(tag,sex) #all tags with their sex

        traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
        tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap

        #VERY LARGE df of all tags for all traps
        #new column of Det.obs - whether that animal was observed in that trap
        #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
        fulldata <- cbind(tags_rep, traps_rep) %>%
          left_join(data, by=c("tag", "trap", "x", "y")) %>%
          group_by(tag) %>% fill(sex, .direction="downup") %>%
          mutate(Det.count = replace_na(Det.count, 0),
                 Det.obs = ifelse(Det.count==0, 0, 1)) %>%
          arrange(tag, trap) %>%
          select(tag, sex, x, y, Det.obs, Det.count) %>%
          rename(x.trap = x,
                 y.trap = y,
                 Tag_ID = tag,
                 Sex = sex)

        #----------------------- 3. Generating overlap network  ------------------------------#

        # Following code from Wanelik & Farine 2022
          # Wanelik, K.M., Farine, D.R. A new method for characterising shared space use networks using animal trapping data.
          # Behav Ecol Sociobiol 76, 127 (2022). https://doi.org/10.1007/s00265-022-03222-5
        
        ##SEASONAL CENTROIDS to calc SEASONAL PARAMS by season, trt, sex, breeder

        # Calculating (weighted) centroids across the SEASON
        # Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
        matrix_dists_real2 <- fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
        matrix_dists_real2 <- matrix_dists_real2 %>% mutate(Det.count2 = Det.count)
        matrix_dists_real2$Det.count2[matrix_dists_real2$Det.count==0]<-1
        matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count2)) #has as many entries (rows) per tag/trap as the animal was captured there
        matrix_dists_real2<-matrix_dists_real2[,-7] #drops the Det.count2 column
        x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
        y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
        obs_centroids <- data.frame(Tag_ID=tagsex$tag, Sex=tagsex$sex)
        obs_centroids$x <- x[match(obs_centroids$Tag_ID,names(x))]
        obs_centroids$y <- y[match(obs_centroids$Tag_ID,names(y))]

        # Creating a matrix of distances between each individual's observed centroid and each trap
        matrix_dists_obs <- fulldata #only one entry per individual
        matrix_dists_obs$x.ind <- obs_centroids$x[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #x coord of centroid of individual
        matrix_dists_obs$y.ind <- obs_centroids$y[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #y coord of centroid of individual
        matrix_dists_obs$Dist <- sqrt((matrix_dists_obs$x.ind-matrix_dists_obs$x.trap)^2 + (matrix_dists_obs$y.ind - matrix_dists_obs$y.trap)^2)
        matrix_dists_obs$Dist.log<-log(matrix_dists_obs$Dist + 1)
        # # Removing unobserved individuals
        # matrix_dists_obs <- matrix_dists_obs[which(!is.na(matrix_dists_obs$Dist)),] #not needed, matrix is only of voles capped that month

        # Calculating total number of detections per individual
        ### TBH Janine isn't sure this is needed because it doesn't appear to be used anywhere, but a good check I guess
        # a <- tapply(matrix_dists_obs$Det.count, list(matrix_dists_obs$Tag_ID),sum)
        # matrix_dists_obs$Det.count.total <- a[match(matrix_dists_obs$Tag_ID,names(a))]

        # Estimating home range profiles (negative sigmoidal curves) using this data
        # Example: Sex-specific curves (male home range profile = fit.males; female home range profile = fit.females)
        fit.males <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Sex=="M"),], family=binomial, control = list(maxit = 50))
        fit.females <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Sex=="F"),], family=binomial, control = list(maxit = 50))

        # Generating overlap network
        params <- data.frame(sex=c("M","F"),
                             a=c(coef(fit.males)[1],coef(fit.females)[1]),
                             b=c(coef(fit.males)[2],coef(fit.females)[2]))

        params <- params %>% mutate(season=season_lab, season_breeder=breeder_lab, trt=trt_lab)

        params_list[[i]][[k]][[j]] <- params

        }

    }

  }

  #collate results
  #this collapses the 2nd order elements (params per trt) down to the 1st order element (season)

  #make a list to store things
  params_list_summary1 <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(params_list)){

    #for both seasons
    summary <- do.call("rbind", params_list[[i]])
    params_list_summary1[[i]] <- summary
  }

  params_list_summary <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(params_list_summary1)){

    #for both seasons
    summary <- do.call("rbind", params_list_summary1[[i]])
    params_list_summary[[i]] <- summary
  }


  #collapse params_list_summary into a df
  params_summary <- do.call(rbind.data.frame, params_list_summary) %>%
    unite("stsb", season, trt, sex, season_breeder)
  row.names(params_summary) <- NULL

  ##### OUTPUT: PARAMS_SUMMARY has the "a" and "b" parameters, calculated per season, per treatment, per sex
  ### params_summary has columns "stsb" (season,treatment,sex,breeding status) "a" and "b"

  #save to RDS file
  saveRDS(params_summary, file = here(params_file))

}

######------------------------------------------------------------------------


## Function for creating spatial overlap networks for individual voles for all the sites, months in a given year
## Inputs: data = FULL fulltrap dataframe that contains all capture data for a given year
##          params_file = file name and extension in "" for the params_file generated by generate_params function
##          centroids_file = file name and extension in "" for the output file of monthly centroids for all voles
##          networks_file = file name and extension in "" for the output file of adjacency matrices
## Output: list (saved to rds file) of all sites, months and an adjacency matrix for each of overlaps between voles
## Output: monthly centroid locations (x,y) for all voles - to be used to create the circle plots by month


create_overlap_networks <- function(data, params_file, centroids_file, networks_file){

  ##---------------- LOAD THE CAPTURE DATA ----------------------

  #load, clean  the fulltrap dataset
  fulltrap <- data %>%
    unite("stsb", season, trt, sex, season_breeder, remove=FALSE) #add stsb column to match params_summary

  #pull all the traps and their x, y coordinates, save as df with (trapID, x, y)
  traps <- fulltrap %>% group_by(trap) %>% slice(1) %>%
    select(trap, x, y) %>%
    arrange(trap)

  #save a vector of the site names (alphabetical order)
  site_names <- unique(fulltrap$site) %>% sort()


  ##---------------- CREATE A LIST of CAPTURE DATA (nested by site, month) --------------------

  #split() makes a list consisting of individual data.frames based on a condition ('site' in this case)
  site_list <- split(fulltrap, fulltrap$site)

  #create new list to hold nested site, month capture
  sitemonth_list <- list()

  for(i in 1:length(site_list)){
    # print(i)
    temp_list <- split(site_list[[i]], site_list[[i]]$month)
    sitemonth_list[[i]] <- temp_list
  }

  #name 1e list element as site (2e list elements are months May-Oct)
  names(sitemonth_list) <- site_names

  ### OUTPUT: SITEMONTH_LIST is a nested list of length 12
  #1e level is all the sites (12)
  #2e level is all the months per site (5) excluding May

  ##-------------- LOAD ADDITIONAL DATA ----------------------

  # Reading in required functions in wanelik_farine_2022_functions.R
  source(here("01-1_wanelik_farine_2022_functions.R"))
  
  # Wanelik, K.M., Farine, D.R. A new method for characterising shared space use networks using animal trapping data.
  # Behav Ecol Sociobiol 76, 127 (2022). https://doi.org/10.1007/s00265-022-03222-5

  #read in the params RDS created in generate_params() function
  params_summary <- readRDS(file = here(params_file))

  ##--------------- (in a loop) GENERATE OVERLAP NETWORK ---------------------

  overlap_network_list <- list()

  centroids_list <- list()

  for(i in 1:length(sitemonth_list)){

    print(names(sitemonth_list[i])) #print the site name

    overlap_network_list[[i]] <- list()

    centroids_list[[i]] <- list()

    for(j in 1:5){
      #1:5 for the five months (MAY DATA OMITTED)

      print(names(sitemonth_list[[i]][j])) #print the month

      data <- sitemonth_list[[i]][[j]] %>%
        select(tag, stsb, trap, x, y) %>%
        group_by(tag, trap) %>%
        mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
        ungroup()

      tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
      tagstsb <- data %>% group_by(tag) %>% slice(1) %>% select(tag, stsb) #all tags with their STSB (season/trt/sex/breeder)

      # skip creating network if there are 0 or 1 animals
      if(nrow(tags)=="0" | nrow(tags)=="1") {next}

      traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
      tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) #all tags replicated for each trap

      #creates VERY LARGE df of all tags for all traps
      #new column of Det.obs - whether that animal was observed in that trap
      #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
      fulldata <- cbind(tags_rep, traps_rep) %>%
        left_join(data, by=c("tag", "trap", "x", "y")) %>%
        group_by(tag) %>% fill(stsb, .direction="downup") %>%
        mutate(Det.count = replace_na(Det.count, 0),
               Det.obs = ifelse(Det.count==0, 0, 1)) %>%
        arrange(tag, trap) %>%
        select(tag, stsb, x, y, Det.obs, Det.count) %>%
        rename(x.trap = x,
               y.trap = y,
               Tag_ID = tag)

      #----------------------- 3. Generating overlap network  ------------------------------#

      # Following code from Wanelik & Farine 2022
        # Wanelik, K.M., Farine, D.R. A new method for characterising shared space use networks using animal trapping data.
        # Behav Ecol Sociobiol 76, 127 (2022). https://doi.org/10.1007/s00265-022-03222-5
      
      ### MONTHLY CENTRIODS with SEASONAL A B PARAMS

      # Recalculating (weighted) centroids - CENTROID OF ALL CAPTURES IN A MONTH
      # Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
      matrix_dists_real2 <- fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
      matrix_dists_real2 <- matrix_dists_real2 %>% mutate(Det.count2 = Det.count)
      matrix_dists_real2$Det.count2[matrix_dists_real2$Det.count==0]<-1
      matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count2)) #has as many entries (rows) per tag/trap as the animal was captured there
      matrix_dists_real2<-matrix_dists_real2[,-7] #drops the Det.count2 column
      x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
      y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
      obs_centroids <- data.frame(Tag_ID=tagstsb$tag, stsb=tagstsb$stsb)
      obs_centroids$x <- x[match(obs_centroids$Tag_ID,names(x))]
      obs_centroids$y <- y[match(obs_centroids$Tag_ID,names(y))]

      # Creating a matrix of distances between each individual's observed centroid and each trap
      matrix_dists_obs <- fulldata #only one entry per individual
      matrix_dists_obs$x.ind <- obs_centroids$x[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #x coord of centroid of individual
      matrix_dists_obs$y.ind <- obs_centroids$y[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #y coord of centroid of individual
      matrix_dists_obs$Dist <- sqrt((matrix_dists_obs$x.ind-matrix_dists_obs$x.trap)^2 + (matrix_dists_obs$y.ind - matrix_dists_obs$y.trap)^2)
      matrix_dists_obs$Dist.log<-log(matrix_dists_obs$Dist + 1)

      # # Removing unobserved individuals
      # matrix_dists_obs <- matrix_dists_obs[which(!is.na(matrix_dists_obs$Dist)),] #not needed, matrix is only of voles capped that month

      # Calculating total number of detections per individual
      # a <- tapply(matrix_dists_obs$Det.count, list(matrix_dists_obs$Tag_ID), sum)
      # matrix_dists_obs$Det.count.total <- a[match(matrix_dists_obs$Tag_ID, names(a))]

      # Estimating home range profiles (negative sigmoidal curves) done separately by season/trt -> params_summary file

      # Generating overlap network
      as <- params_summary$a[match(obs_centroids$stsb,params_summary$stsb)]
      bs <- params_summary$b[match(obs_centroids$stsb,params_summary$stsb)]

      overlap_network_sep <- get_network_2D(obs_centroids$x[which(!is.na(obs_centroids$x))], obs_centroids$y[which(!is.na(obs_centroids$y))],as,bs)
      rownames(overlap_network_sep)=colnames(overlap_network_sep)=na.omit(obs_centroids)$Tag_ID

      overlap_network_list[[i]][[j]] <- overlap_network_sep

      #save monthly centroids (for plotting)
      centroids_list[[i]][[j]] <- obs_centroids %>% select(Tag_ID, x, y) %>% mutate(month=paste(names(sitemonth_list[[i]][j])),
                                                                                    site=paste(names(sitemonth_list[i])))

    }


  }


  #name the 12 1st order elements of overlap_network_list as the sites
  names(overlap_network_list) <- names(sitemonth_list)

  #rename the sublist items (months) for each site
  for(i in 1:length(overlap_network_list)){
    names(overlap_network_list[[i]]) <- c("june", "july", "aug", "sept", "oct")
  }

  #in case there were NULL list elements (for a month with 0 or 1 animals)
  library(rlist)
  #a slick little function that removes NULL list elements (recursive=TRUE) to work through nested lists
  overlap_network_list <- list.clean(overlap_network_list, fun = is.null, recursive = TRUE)

  #save networks file since that code takes 5ever to run
  saveRDS(overlap_network_list, here(networks_file))

  

  #collate MONTHLY centroids results
  #make a list to store things
  centroids_summary1 <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(centroids_list)){

    #for all 12 sites
    summary <- do.call("rbind", centroids_list[[i]])
    centroids_summary1[[i]] <- summary
  }

  #name the 12 1st order elements as their sites
  names(centroids_summary1) <- names(sitemonth_list)

  ## make net_mets_list_summary into freiggein huge df
  centroids_summary <- do.call(rbind.data.frame, centroids_summary1)
  row.names(centroids_summary) <- NULL

  #saveRDS
  saveRDS(centroids_summary, here(centroids_file))

}

#####-------------------------------------------------------------------------------------

#### THIS CODE calculates the network metrics (mostly just weighted degree and its relatives)
# for all the voles at all the times
# output is "wt_net_mets_summary" file which can be used for downstream analysis

## Input: data = the FULL fulltrap df for that year
##        networks_file = networks file generated in create_overlap_networks
##        netmets_file = file to be generated, df of network metrics
## Output: netmets_file = dataframe of network metrics for every vole in every occasion it was captured


calculate_network_metrics <- function(data, networks_file, netmets_file){

  ##---------------- LOAD THE DATA ----------------------

  #load, clean the fulltrap dataset
  fulltrap <- data %>%
    unite("sts", season, trt, sex, remove=FALSE) #add sts column to match params_summary

  # load the network data
  overlap_network_list <- readRDS(here(networks_file))


  ##--------------------------------------------

  # Create list to store results
  wt_net_mets_list <- list()

  # Calculate network metrics
  for(i in 1:length(overlap_network_list)){
    #using 'length(overlap_network_list)' to loop through each of the 12 sites

    #for each site
    print(names(overlap_network_list[i]))
    site.id <- names(overlap_network_list[i])

    site <- list()

    for(j in 1:length(overlap_network_list[[i]])){
      #using 'length(overlap_network_list[[i]]))' to loop through all trapping occasions (ie months) with capture data per site

      #for each month
      print(names(overlap_network_list[[i]][j]))
      month.id <- names(overlap_network_list[[i]][j])

      adjmat <- overlap_network_list[[i]][[j]]

      #create WEIGHTED NETWORK from adjacency matrix
      inet <- graph_from_adjacency_matrix(adjmat, mode="undirected", weighted = TRUE, diag = FALSE)

      #create UNWEIGHTED NETWORK (just for binary degree)
      adjmatscaled <-ifelse(adjmat>0.01,1,0)
      inet_bin01 <- graph_from_adjacency_matrix(adjmatscaled, weighted=NULL, mode="undirected", diag=FALSE)

      #and for SUPPLEMENT, a lil sensitivity testing
      adjmatscaled <-ifelse(adjmat>0.05,1,0)
      inet_bin05 <- graph_from_adjacency_matrix(adjmatscaled, weighted=NULL, mode="undirected", diag=FALSE)

      adjmatscaled <-ifelse(adjmat>0.001,1,0)
      inet_bin001 <- graph_from_adjacency_matrix(adjmatscaled, weighted=NULL, mode="undirected", diag=FALSE)

      adjmatscaled <-ifelse(adjmat>0.005,1,0)
      inet_bin005 <- graph_from_adjacency_matrix(adjmatscaled, weighted=NULL, mode="undirected", diag=FALSE)
      ##end sensitivity testing##

      ids <- get.vertex.attribute(inet, "name") #tag ids for all the animals on the grid
      month <- rep(names(overlap_network_list[[i]])[j],length(ids)) #capture month

      #dataframe to hold results per month
      site[[j]] <- data.frame(ids, month)

      #network metrics to calculate
      site[[j]]$wt.deg <- igraph::strength(inet) #this is the sum of all degree weights for a node
      site[[j]]$norm.wt.deg <- (igraph::strength(inet))/((igraph::gorder(inet))-1) #your strength/(total nodes-you)

      #binary degree (number of overlaps)
      site[[j]]$bin.01.deg <- igraph::degree(inet_bin01)

      site[[j]]$bin.05.deg <- igraph::degree(inet_bin05)
      site[[j]]$bin.001.deg <- igraph::degree(inet_bin001)
      site[[j]]$bin.005.deg <- igraph::degree(inet_bin005)

      #number of nodes (voles) in the network
      site[[j]]$n.node <- rep(igraph::gorder(inet), length(ids))

    }

    #write the list 'site' as a 1st-order item in wt_net_mets_list
    wt_net_mets_list[[i]] <- site

  }

  #name the 12 1st order elements of wt_nets_list as the sites
  names(wt_net_mets_list) <- names(overlap_network_list)

  #rename the sublist items (months) for each site
  for(i in 1:length(wt_net_mets_list)){
            names(wt_net_mets_list[[i]]) <- c("june", "july", "aug", "sept", "oct")
  }


  ################################### ABOUT wt_net_mets_list #############################################
  #the output of wt_net_mets_list is a list of 12 1st-order items, 1 per site
  #under each site, there are 5 second-order items, 1 per trapping occasion (named by month)
  #each of those 2nd-order items (months) is a df with tag ID, month, and all the network metrics...
  #for ONLY THE ANIMALS captured on that grid, during that occasion
  ########################################################################################################


  ############## condense wt_net_mets_list to make it easier to use for analysis #######################

  #collate results
  #this collapses the 2nd order elements (network metrics for a single month) down to the 1st order element (site)
  #so now wt_net_mets_list_summary is a list of 12 dfs, each is all the network metrics for a site across all the months

  #make a list to store things
  wt_net_mets_list_summary <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(wt_net_mets_list)){

    #for all 12 sites
    summary <- do.call("rbind", wt_net_mets_list[[i]])
    wt_net_mets_list_summary[[i]] <- summary
  }

  #remove the row names of the dfs (its the month, but we already have that info)
  for(i in 1:length(wt_net_mets_list_summary)){
    row.names(wt_net_mets_list_summary[[i]]) <- NULL
  }

  #name the 12 1st order elements as their sites
  names(wt_net_mets_list_summary) <- names(overlap_network_list)

  ## make net_mets_list_summary into freiggein huge df
  wt_net_mets_summary <- do.call(rbind.data.frame, wt_net_mets_list_summary)

  #clean up the df
  wt_net_mets_summary <- wt_net_mets_summary %>%
    rownames_to_column("name") %>% #row names are the sites, make that a column
    separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
    mutate(site = as.factor(site)) %>% #make site a factor
    rename(tag=ids)

  #save it
  saveRDS(wt_net_mets_summary, here(netmets_file))


}


