# Pseudotime estimation with Monocole2 - ROS/MAP brain metabolomics monocle2 analysis


setwd('~/multiomic_pseudotimes/')


# clear workspace
rm(list=ls())

# libraries
BiocManager::install("SummarizedExperiment")
devtools::install_github(repo="krumsieklab/maplet@v1.1.2", subdir="maplet")
library(maplet) # omics analysis toolbox
library(tidyverse) # tidy
library(magrittr) # tidy %<>%
library(monocle) # pseudotime
library(readxl) # multiple excel sheets
library(openxlsx) # excel
#source('/Users/rib4003/Projects/rosmap/pseudotime/custom_functions.R')

# input - preprocessed brain metabolomics dataset
#rosmap_met_pdata <- 'data/tmp_rosmap_brain_metabolomics_processed_medcor_data.xlsx'

metsFid <- "syn48321418"
synapser::synGet(metsFid, downloadLocation = "files/")

path <- readxl_example("files/rosmap_metabolomics_processed_med_covar_corrected_females.xlsx")

path <- "files/rosmap_metabolomics_processed_med_covar_corrected_females.xlsx"
sheetnames <- excel_sheets(path)

mylist <- lapply(excel_sheets(path), read_excel, path = path)

names(mylist) <- sheetnames

assay <- mylist$assay
metabolite_metadata <- mylist$rowData
metadata <- mylist$colData

# load preprocessed medication corrected data
Dmc <- mt_load_se_xls(file=metF$path) %>%
  #inverse cognitive decline slope for consistency
  mt_anno_mutate(anno_type = "samples", col_name = 'cogng_random_slope', 
                 term = (-1*cogng_random_slope))

# turn outcomes and covariates into appropriate data types factor / numeric
Dmc %<>% mt_anno_mutate(anno_type = "samples", col_name = 'msex', 
                        term = as.factor(msex)) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'apoe_genotype', 
                 term = as.factor(apoe_genotype)) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'pmi', 
                 term = as.numeric(as.matrix(pmi))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'age_death', 
                 term = as.numeric(as.matrix(age_death))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'bmi', 
                 term = as.numeric(as.matrix(bmi))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'educ', 
                 term = as.numeric(as.matrix(educ))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'niareagansc', 
                 term = as.factor(niareagansc)) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'sqrt_amyloid', 
                 term = as.numeric(as.matrix(sqrt_amyloid))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'tangles', 
                 term = as.numeric(as.matrix(tangles))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'cogng_random_slope', 
                 term = as.numeric(as.matrix(cogng_random_slope))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'cogn_global', 
                 term = as.numeric(as.matrix(cogn_global))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'gpath', 
                 term = as.numeric(as.matrix(gpath))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'cogdx', 
                 term = as.factor(as.matrix(cogdx))) %>%
  # sage's diagnosis definition
  # Control: cogdx == 1, braaksc <= 3, and ceradsc >= 3
  #AD: cogdx == 4, braaksc >= 4, and ceradsc <= 2
  mt_anno_mutate(anno_type = "samples", col_name = 'sage_diag', 
                 term = case_when((cogdx==1 & diagnosis==0) ~0,  
                 (cogdx==4 & diagnosis==1) ~1, TRUE ~ NA_real_)) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'sage_diag', 
                 term = as.factor(as.matrix(sage_diag))) %>%
  {.}

# correct for covariates
D1 <- mt_pre_confounding_correction(Dmc, formula = '~ age_death + msex + pmi + bmi + educ + apoe_genotype')

#upload sex-specific matrices
F1 
# association analysis with sage's diagnosis
D3 <- mt_stats_univ_lm(D1, formula = ~sage_diag, stat_name = 'sage_diag') %>% 
  mt_post_multtest(stat_name = 'sage_diag', method='BH')
S <- D3 %>% maplet:::mtm_res_get_entries("stats")
stat_res <- S[[1]]$output$table %>% data.frame()

# monocole2 estimation all samples ----
# set seed for reproducibility
set.seed(484956)
## monocle data object with all metabolites
## unsupervised
hsmm_usup <- get_trajectories(D1)
## monocle data object with all metabolites significantly associated with sage's diagnosis
D2 <- D1 %>% mt_modify_filter_features(filter=make.names(rowData(D1)$name, unique = T)%in%(stat_res %>% filter(p.adj<0.05) %>% pull(var)))
## supervised
hsmm_sup <-  get_trajectories(D2)
# check to see if State designations make sense
## unsupervised
p_usup <- plot_cell_trajectory(hsmm_usup, color_by = "State",
                          show_branch_points=F,
                          use_color_gradient = F,cell_size = 1.5) + 
  ggplot2::scale_color_viridis_d() +
  ggplot2::labs(color="State")
## supervised
p_sup <- plot_cell_trajectory(hsmm_sup, color_by = "State",
                          show_branch_points=F,
                          use_color_gradient = F,cell_size = 1.5) + 
  ggplot2::scale_color_viridis_d() +
  ggplot2::labs(color="State")
###reverse supervised trajectory
hsmm_sup <- orderCells(hsmm_sup, reverse = TRUE)

# outcomes to project on the trajectories
outcomes <- read.xlsx('/Users/rib4003/Projects/rosmap/pseudotime/data/outcomes.xlsx')
## loop over outcomes for unsupervised plots
plot_usup <- get_plot_for_all_outcomes(hsmm=hsmm_usup, outcomes=outcomes)
save(plot_usup, file='output_usup.rds')
list_to_rmd(plot_usup, outfile = "output_usup.html")

## loop over outcomes for supervised plots
plot_sup <- get_plot_for_all_outcomes(hsmm=hsmm_sup, outcomes=outcomes)
save(plot_sup, file='output_sup.rds')
list_to_rmd(plot_sup, outfile = "output_sup.html")
