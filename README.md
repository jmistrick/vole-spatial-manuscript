# vole-spatial-manuscript
Associated code and data files for manuscript: "Effects of Food Supplementation and Helminth Removal on Space Use and Spatial Overlap in Wild Rodent Populations"

J. Mistrick, J.S.M. Veitch, S.M. Kitchen, S. Clague, B.C. Newman, R.J. Hall, S.A. Budischak, K.M. Forbes, M.E. Craft

In this study, we leverage a large-scale field experiment of wild bank voles (Clethrionomys glareolus) consisting of resource supplementation and parasite removal by oral dewormer to quantify the effects of food supplementation and helminth removal on space use and spatial overlap in wild rodent populations.

A series of R scripts and data files are used to estimate space use following the methods of Wanelik & Farine (2022) and construct spatial overlap networks in the vole populations. The analysis conducted here and figures and tables generated all appear in the aforementioned manuscript.

The files are numbered to direct the reader to their sequence for analysis:
1. Analysis of factors (sex, reproductive status, treatment) influencing seasonal space use (01_factors_influencing_seasonal_space_use.R). This script sources the functions developed by Wanelik & Farine (2022) which can be found in the file (01-1_wanelik_farine_2022_functions.R).
2. Estimating bank vole space use, constructing spatial overlap networks, quantifying individual spatial overlap between pairs of bank voles (02_construct_spatial_overlap_networks.R). This script sources functions developed by Wanelik & Farine (2022) and modified for this analysis by Janine Mistrick found in the file (02-1_functions_construct_overlap_networks).
3. Data analysis and visualization including all figures and tables found in the main text and supplemental materials (03_analysis_tables_figures.R).
4. Brief analysis of the efficacy of the deworming treatment (04_deworm_efficacy.R)


Wanelik, K.M., Farine, D.R. A new method for characterising shared space use networks using animal trapping data.Behav Ecol Sociobiol 76, 127 (2022). https://doi.org/10.1007/s00265-022-03222-5
