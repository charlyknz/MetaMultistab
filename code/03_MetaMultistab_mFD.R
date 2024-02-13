#### R script for meta analysis ####

rm(list=ls())
graphics.off()

library(tidyverse)
library(readxl)
library(here)
library(cowplot)
library(GGally)
library(ggpubr)

#### import data ####

# species stability
SpeciesStab <- read_csv('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/SpeciesStabilities.csv')%>%  select(-'...1')

# community stability
ComStab <- read_csv('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/CommunityStabilities.csv')%>%  select(-'...1')

# Merged Stability
MergedComStab<- read_csv('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/MergedCommunityStability.csv') %>%  select(-'...1')


#### mFD - traits ####
library(mFD)


## Data requirements: NAs are 0! 
## multiple assemblages to compare among 
## occurrences or abundance data to tell which species are in one assemblage (minimum presence absence)
## sp_tr = trait data (AUC.pi and RR) needs: species as column names 
## asb_sp_w = community matrix based on presence absence data in control. This is needed to see which species are co-occurring
## as species need to be unique: here we use populatin level instead <- combination of caseID and species


### Here we have continuous traits only, therefore follow: 
# https://cmlmagneville.github.io/mFD/articles/Continuous_traits_framework.html ###

#' Compute Functional Identity
#'
#' This function computes the weighted average position along each axis. FIde 
#' is computed using relative weight so that it is not affected by unit 
#' (e.g. g or kg for biomass). In the special case where 'weight' is filled 
#' with only 0/1 (absence/presence), then FIde will be computed assuming that
#' all species have the same weight. The results of this function are used in 
#' FSpe, FOri and FNND computation.
#'
#' @param sp_faxes_coord_k a matrix of species coordinates present in a
#'   given assemblage in a chosen functional space with only needed axes.
#'   Species coordinates have been retrieved thanks to \code{tr.cont.fspace} or
#'   \code{\link{quality.fspaces}} and filtered thanks to
#'   \code{\link{sp.filter}}.
#'
#' @param asb_sp_relatw_k a matrix containing species relative weight
#'   (columns) for a given assemblage.
#'   
#' @param k a character string referring to the assemblage studied.
#'
#' @param check_input a logical value allowing to test or not the
#'   inputs. Possible error messages will thus may be more understandable for
#'   the user than R error messages. Species coordinates matrix and
#'   species*weight data frame must not contain NA, their rownames must be 
#'   filled and they must have similar names values. 
#'   Default: `check_input = FALSE`.
#'
#' @return A matrix containing functional identity values for a given 
#'   assemblage along the dimensions (columns). Number of dimensions is fixed 
#'   to the number of dimensions in \code{sp_faxes_coord} data frame.
#'   
data.plot <-SpeciesStab

# Trait data
TraitData <- data.plot %>%
  group_by(caseID) %>%
  mutate(Con.M = as.numeric((sum.con))) %>%
  filter(Con.M >0)
TraitData$species<-gsub( ' ', '_',TraitData$species)

sp_tr <- TraitData %>%
  ungroup()%>%
  filter(!caseID %in%  c('CK015_5','CK015_6','CK015_7','CK015_8','CK015_9','CK015_15',"CK015_10",'CK015_20','CK015_21','CK015_22','CK015_23',
                         'CK027_1','CK027_3','CK041_1','CK041_10','CK041_11','CK041_12','CK041_2','CK041_3','CK041_4','CK041_5',
                         'CK041_6','CK041_7','CK041_8','CK041_9',
                         'HH002_1')) %>%
  mutate(population = paste(species, caseID, sep = '_')) %>%
  select(population, AUC.pi, AUC.RR)%>%
  column_to_rownames(var = 'population') 
sp_tr[is.na(sp_tr)] <-0


# community matrix
asb_sp_w <- TraitData %>%
  filter(!caseID %in%  c('CK015_15','CK015_5','CK015_6','CK015_7','CK015_8','CK015_9',"CK015_10",'CK015_20','CK015_21','CK015_22','CK015_23',
                         'CK027_1','CK027_3','CK041_1','CK041_10','CK041_11','CK041_12','CK041_2','CK041_3','CK041_4','CK041_5',
                         'CK041_6','CK041_7','CK041_8','CK041_9',
                         'HH002_1')) %>%
  mutate(population = paste(species, caseID, sep = '_')) %>%
  select(population, caseID, Con.M) %>%
  spread(key = population, value = Con.M) %>%
  column_to_rownames(var = 'caseID') 
asb_sp_w[is.na(asb_sp_w)] <-0
asb_sp_w[(asb_sp_w>0)] <-1
str(asb_sp_w)
asb_sp_w<-as.matrix(asb_sp_w)

# 2. compute functional space
fspace <-mFD::tr.cont.fspace(
  sp_tr        = sp_tr, 
  pca          = FALSE, 
  nb_dim       = 2, 
  scaling      = "scale_center",
  compute_corr = "none")

fspace
summary(fspace)
fspace$"sp_faxes_coord"


# Plot functional space
big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = fspace$"sp_faxes_coord",
  faxes           = c('AUC.RR','AUC.pi'),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

big_plot$patchwork

alpha_fd_indices <- mFD::alpha.fd.multidim( #add  mFD to standardise by global 
  sp_faxes_coord  = fspace$"sp_faxes_coord",
  asb_sp_w         = asb_sp_w,
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric",  
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)


# export functional diversity indices
fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"
fd_ind_values

# look at details
details_list <- alpha_fd_indices$"details"


## plot exemplary two caseID hypervolumina within the maximum volumina 
plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("HH018_1",'CK002_1'),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_vert               = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_sp                  = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_vert                = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_ch                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_ch                  = c(pool = "white", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

plots_alpha$"fric"$"patchwork"

FD_ind_values <- fd_ind_values %>%
  tibble::rownames_to_column() %>%
  rename(caseID = rowname) 

#### Merge Fdiv and Functional community stability ####

#new df
Fdiv <- FD_ind_values%>%
  distinct(caseID,sp_richn,fdis,fmpd,fnnd,feve,fric,fspe)%>%
  left_join(., ComStab, by = c('caseID')) 
which(is.na(Fdiv$AUC.sum.delatbm.tot))
#

#remove
#Fdiv$organism[Fdiv$organism == 'Cladoceran'] <- 'zooplankton'
#Fdiv$organism[Fdiv$organism == 'invertebrate'] <- 'invertebrates'
#Fdiv$organism[Fdiv$organism == 'macroinvertebrate'] <- 'invertebrates'

unique(Fdiv$organism)
#Fdiv$OrganismType <- Fdiv$organism
#Fdiv$OrganismType[Fdiv$OrganismType == 'macroalgae'] <- 'primaryproducer'
#Fdiv$OrganismType[Fdiv$OrganismType == 'phytoplankton'] <- 'primaryproducer'

#Fdiv$OrganismType[Fdiv$OrganismType == 'periphyton'] <- 'primaryproducer'
#Fdiv$OrganismType[Fdiv$OrganismType == 'macrophytes'] <- 'primaryproducer'
#Fdiv$OrganismType[Fdiv$OrganismType == 'plant'] <- 'primaryproducer'
#Fdiv$OrganismType[Fdiv$OrganismType == 'zooplankton'] <- 'herbivore'
#Fdiv$OrganismType[Fdiv$OrganismType == 'vertebrates'] <- 'predator'
#Fdiv$OrganismType[Fdiv$OrganismType == 'fish'] <- 'predator'
#Fdiv$OrganismType[Fdiv$OrganismType == 'invertebrates'] <- 'predator'
#Fdiv$OrganismType[Fdiv$OrganismType == 'invertebrate'] <- 'predator'
#Fdiv$OrganismType[Fdiv$OrganismType == 'macroinvertebrate'] <- 'predator'


str(Fdiv)
#write.csv(Fdiv, file = 'Output/StabFdiv_n4.csv')

### Plot FDiv Indices ~ Stability ###

Fdiv %>%
 # select(caseID, fdis, fric, feve, sp_richn,studyID, organism, dist.cat, OEV, resp.cat,system,OrganismType) %>%
  gather(c(fdis, fric, feve), key = 'Indices', value = 'IndexValue') %>%
  ggplot(., aes(x = IndexValue, y = AUC.sum.delatbm.tot))+
  geom_point(size = 2, alpha = .6)+
  labs(x = 'Fdiv Value', y = 'Community Instability')+
  facet_grid(~Indices, scales = 'free_x') +
  theme_bw()+
  theme(text = element_text(size=rel(4)),
        strip.text.x = element_text(size=rel(4)))
#ggsave(last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/Fdiv_overview.png'), width = 8, height = 5)


sp_richn<-ggplot(Fdiv, aes(x = sp_richn, y = AUC.sum.delatbm.tot))+
  geom_hline(yintercept = 0)+
  geom_point(size = 2, alpha = .4)+
  labs(x = 'S Richness', y = 'AUC.sum.delatbm.tot')+
  theme_bw()+
  theme(legend.position = 'bottom')
sp_richn

ggsave(last_plot(), file = here('~/Desktop/phD/Meta_Multistab/FDiv_LinFdivAdj.png'), width = 14, height = 10)

#### START Meta-Analysis ####

#libraries we need
library(metafor)
library(broom)

## look at Effect sizes 

study <- read_excel("~/Desktop/phD/Meta_Multistab/MetaMultistab/Multistab_species_data_mfd.xlsx") %>%
  select(-c(35:51)) %>%
  rename(caseID = paper)

metadata_Fdiv <- left_join(Fdiv, study, by = c('caseID')) %>%
  distinct(caseID ,studyID,duration, open, organism, system,dist.cat,OEV,AUC.sum.delatbm.tot, fdis, fric, feve)

hist(metadata_Fdiv$OEV)
unique(metadata_Fdiv$caseID)


# unweighted MA requires column 1
metadata_Fdiv$unweighted<-1
names(metadata_Fdiv)

Fdiv.m<-rma.mv(AUC.sum.delatbm.tot,unweighted,
                      mods = ~fdis+fric+feve+dist.cat+system+open+organism+duration,
                      random = ~ 1 | caseID,
                      method="REML",data=metadata_Fdiv)
summary(Fdiv.m) #175.1314
complete.mod<-tidy(summary(Fdiv.m))
complete.mod$k<-Fdiv.m$k
complete.mod$AIC<-Fdiv.m$fit.stats$REML[3]

#2
Fdiv.m<-rma.mv(AUC.sum.delatbm.tot,unweighted,
               mods = ~fdis+fric+feve+system,
               random = ~ 1 | caseID,
               method="REML",data=metadata_Fdiv)
summary(Fdiv.m) #174.2335 
