#### R script for meta analysis ####

rm(list=ls())
graphics.off()

library(tidyverse)
library(readxl)
library(here)
library(cowplot)
library(GGally)
library(ggpubr)
library(ggpmisc)


#### import data ####

# species stability
SpeciesStab <- read_csv('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/SpeciesStabilities.csv')%>%  select(-'...1')

# community stability
ComStab <- read_csv('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/CommunityStabilities.csv')%>%  select(-'...1') %>%
  rename(Deltabm.Tot = AUC.sum.delatbm.tot,
         Recovery=Recov,
         Resistance=Resist)

# Merged Stability
MergedComStab<- read_csv('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/MergedCommunityStability.csv') %>%  select(-'...1')
names(MergedComStab)


#### Correlation community stability ####
ggscatter(MergedComStab,  x= 'AUC.delatbm.tot.MA', y = 'AUC.sum.delatbm.tot', add = 'reg.line',cor.coef = TRUE)
ggsave(plot= last_plot(), file = 'Correlation_MA_sum_Stability.png')

#### ResponseDiversity ####

source(here("~/Desktop/phD/Meta_Multistab/response-diversity-pulse-pert/R/0-functions/Ross_et_al_functions.R"))

names(SpeciesStab)

realised.pert <- SpeciesStab %>%
  group_by(caseID, resp.cat) %>%
  reframe(mean_spp_deltabm = mean(AUC.RR),
            var_spp_deltabm = var(AUC.RR),
            RD_diss = resp_div(AUC.RR, sign_sens = FALSE),
            RD_div = resp_div(AUC.RR, sign_sens = TRUE),
            mean_abs_spp_deltabm = mean(abs(AUC.RR)),
            var_abs_spp_deltabm = var(abs(AUC.RR)),
            RD_abs_diss = resp_div(abs(AUC.RR), sign_sens = FALSE),
            RD_abs_div = resp_div(abs(AUC.RR), sign_sens = TRUE)) %>%
  ungroup()%>%
  gather(mean_spp_deltabm,var_spp_deltabm,RD_diss,RD_div, mean_abs_spp_deltabm,var_abs_spp_deltabm,RD_abs_div,RD_abs_diss,key = 'RD.metric', value = 'RD.value') 

#### Community stability and RD ####

##### Sum CommunityStability ####
AllStab <- merge(realised.pert,ComStab, by = c('caseID', 'resp.cat')) %>%
  filter(!str_detect(RD.metric,'abs' ))

str(AllStab)


AllStab_wide <-  AllStab %>%
  spread(key = RD.metric, value = RD.value)  
mean1<- ggscatter(AllStab_wide, x = 'mean_spp_deltabm', y = 'Deltabm.Tot', 
                  xlab = 'Mean species response trait', ylab = 'Relative OEV',
                  add = 'reg.line', cor.coef = T, cor.method='spearman')
mean1

RD_div1<- ggscatter(AllStab_wide, x = 'RD_div', y = 'Deltabm.Tot', 
                  xlab = 'Response divergence', ylab = 'Relative OEV',
                  add = 'reg.line', cor.coef = T, cor.method='spearman')
RD_div1

RD_diss1<- ggscatter(AllStab_wide, x = 'RD_diss', y = 'Deltabm.Tot', 
                    xlab = 'Response dissimilarity', ylab = 'Relative OEV',
                    add = 'reg.line', cor.coef = T, cor.method='spearman')
RD_diss1

plot_grid(RD_div1, RD_diss1,mean1, labels = c('(a)', '(b)', '(c)'))
ggsave(plot = last_plot(), file = here('output/Correlation_RealisedResponseTraits_Instab.png'), width = 8, height = 7)


### Stability Metrics ###
AllStab %>%
  filter(!str_detect(RD.metric,'abs' ))%>%
  ggplot(., aes ( x = RD.value, y = Resistance))+
  geom_point()+
  # geom_smooth(method = lm, se =F, alpha = 0.6)+
  facet_wrap(~RD.metric, scales = 'free_x')+
  theme_bw()+
  labs(x = 'Realised Response Traits', y = 'Resistance')
ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/ResponseTraits_Resist_merged.png'), width =  6, height = 5)

AllStab %>%
  filter(!str_detect(RD.metric,'abs' ))%>%
  ggplot(., aes ( x = RD.value, y = Recovery))+
  geom_point()+
  # geom_smooth(method = lm, se =F, alpha = 0.6)+
  facet_wrap(~RD.metric, scales = 'free_x')+
  theme_bw()+
  labs(x = 'Realised Response Traits', y = 'Recovery')
ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/ResponseTraits_Recov_merged.png'), width =  6, height = 5)


AllStab %>%
  filter(!str_detect(RD.metric,'abs' ))%>%
  ggplot(., aes ( x = RD.value, y = CV))+
  geom_point()+
  # geom_smooth(method = lm, se =F, alpha = 0.6)+
  facet_wrap(~RD.metric, scales = 'free_x')+
  theme_bw()+
  labs(x = 'Realised Response Traits', y = 'Temporal Variability (CV)')
ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/ResponseTraits_CV_merged.png'), width =  6, height = 5)


##### Merge RD and community stability from sum #####
str(ComStab)
AllStab1 <- merge(igr.pert,ComStab, by = c('caseID', 'resp.cat'))
str(AllStab1)

p1 <- AllStab1 %>%
  filter(!str_detect(RD.metric,'abs' ))%>%
  ggplot(., aes ( x = RD.value, y = Deltabm.Tot))+
  geom_point()+
  labs(x = 'Realised Response Traits', y = 'rel. OEV')+
  facet_wrap(~RD.metric, scales = 'free_x')+
  theme_bw()
p1
ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/ResponseTraits_relOEV_sum.png'), width =  6, height = 5)


#### Stability Metrics ####
p2<-AllStab1 %>%
  filter(!str_detect(RD.metric,'abs' ))%>%
  ggplot(., aes ( x = RD.value, y = Resistance))+
  geom_point()+
  # geom_smooth(method = lm, se =F, alpha = 0.6)+
  facet_wrap(~RD.metric, scales = 'free_x')+
  theme_bw()+
  labs(x = 'Realised Response Traits', y = 'Resistance')
p2
ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/ResponseTraits_Resist_sum.png'), width =  6, height = 5)

p3<-AllStab1 %>%
  filter(!str_detect(RD.metric,'abs' ))%>%
  ggplot(., aes ( x = RD.value, y = Recovery))+
  geom_point()+
  # geom_smooth(method = lm, se =F, alpha = 0.6)+
  facet_wrap(~RD.metric, scales = 'free_x')+
  theme_bw()+
  labs(x = 'Realised Response Traits', y = 'Recovery')
p3
ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/ResponseTraits_Recov_sum.png'), width =  6, height = 5)


p4<-AllStab1 %>%
  filter(!str_detect(RD.metric,'abs' ))%>%
  ggplot(., aes ( x = RD.value, y = CV))+
  geom_point()+
  # geom_smooth(method = lm, se =F, alpha = 0.6)+
  facet_wrap(~RD.metric, scales = 'free_x')+
  theme_bw()+
  labs(x = 'Realised Response Traits', y = 'Temporal Variability (CV)')
p4
ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/ResponseTraits_CV_sum.png'), width =  6, height = 5)

plot_grid(p1, p2,p3,p4, labels = c('a)', 'b)', 'c)', 'd)'))
ggsave(plot = last_plot(), file= here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/SumStabilityMetric.png'), width = 10, height = 10)

#### START Meta-Analysis ####

#libraries we need
library(metafor)
library(broom)

## look at Effect sizes 
summary(AllStab)
names(AllStab)
study <- read_excel("~/Desktop/phD/Meta_Multistab/MetaMultistab/Multistab_species_data_mfd.xlsx") %>%
  select(-c(35:51)) %>%
  rename(caseID = paper)

metadata <- left_join(AllStab, study, by = c('caseID', 'resp.cat')) %>%
  drop_na(RD.metric) %>%
  distinct(caseID ,studyID,duration, open, organism, system,dist.cat,resp.cat,OEV,Deltabm.Tot, RD.metric, RD.value) %>%
  filter(!str_detect(RD.metric, 'abs'))%>%
  spread(key = RD.metric, value = RD.value) 

setdiff(AllStab$caseID,metadata$caseID)
hist(metadata$OEV)
unique(metadata$caseID)

metadata <- filter(metadata, resp.cat !=  "contribution to production")

# unweighted MA requires column 1
metadata$unweighted<-1
names(metadata)

#### Delatbm.tot ####
#m0
m0<-rma.mv(Deltabm.Tot,unweighted,
                      mods = ~RD_diss+RD_div+var_spp_deltabm+mean_spp_deltabm+resp.cat,
                      random = ~ 1 | caseID,
                      method="REML",data=metadata)
summary(m0) # 188.3392

## m1
m01<-rma.mv(Deltabm.Tot,unweighted,
                      mods = ~RD_diss+RD_div+var_spp_deltabm+mean_spp_deltabm,
                      random = ~ 1 | caseID,
                      method="REML",data=metadata)
summary(m01)#187.7341

##m2
#remove organism, duration, resp.cat, dist.cat,open,system
m02<-rma.mv(Deltabm.Tot,unweighted,
                      mods = ~mean_spp_deltabm,
                      random = ~ 1 | caseID,
                      method="REML",data=metadata)
summary(m02)#185.8354


#### SUM Deltabm.tot ####
metadata1 <- left_join(AllStab1, study, by = c('caseID', 'resp.cat')) %>%
  drop_na(RD.metric) %>%
  distinct(caseID ,studyID,duration, open, organism, system,dist.cat,resp.cat,OEV,Deltabm.Tot, RD.metric, RD.value) %>%
  filter(!str_detect(RD.metric, 'abs'))%>%
  spread(key = RD.metric, value = RD.value) 

setdiff(AllStab1$caseID,metadata$caseID)
unique(metadata1$caseID)

metadata1 <- filter(metadata1, resp.cat !=  "contribution to production")

# unweighted MA requires column 1
metadata1$unweighted<-1
names(metadata1)

test.complete<-rma.mv(Deltabm.Tot,unweighted,
                      mods = ~RD_diss+RD_div+var_spp_deltabm+mean_spp_deltabm+resp.cat,
                      random = ~ 1 | caseID,
                      method="REML",data=metadata1)
summary(test.complete) #185.9527
complete.mod<-tidy(summary(test.complete))
complete.mod$k<-test.complete$k
complete.mod$AIC<-test.complete$fit.stats$REML[3]

#2
test.complete<-rma.mv(Deltabm.Tot,unweighted,
                      mods = ~var_spp_deltabm+mean_spp_deltabm,
                      random = ~ 1 | caseID,
                      method="REML",data=metadata1)
summary(test.complete)#185.6811
complete.mod<-tidy(summary(test.complete))
complete.mod$k<-test.complete$k
complete.mod$AIC<-test.complete$fit.stats$REML[3]

### forest ###
forest_plot <- forest(test.complete$b, # This is our vector of two effect sizes for our 2 sub-categories
                      ci.lb = test.complete$ci.lb, # Confidence interval lower bounds
                      ci.ub = test.complete$ci.ub, # Confidence interval upper bounds
                      annotate = TRUE, # This tells function to list effect size and CIs for each group on our graph
                      xlab = "ln(Response Ratio)", # label for x-axis
                      slab = c('intrcpt','var_spp_deltabm','mean_spp_deltabm'), #label for y-axis
                      cex = 1, # Font size for entire graph (excluding headers)
                      digits = 2 # Round effect size and CI to 2 digits
)
op <- par(cex=1, font=2) # Set up font for rest of graph (just the headers of the graph remain), to make bold headings, set font=2
text(2.7, 3.5, "ln(Response Ratio) [95% CI]")
text(-2.1,3.5, "Mods")


