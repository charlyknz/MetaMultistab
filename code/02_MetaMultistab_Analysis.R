#### R script for meta analysis ####
graphics.off()

library(tidyverse)
library(readxl)
library(here)
library(cowplot)
library(GGally)
library(ggpubr)
library(ggpmisc)
library(metafor)


#### import data ####

# species stability
SpeciesStab <- read_csv('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/SpeciesStabilities.csv')%>%  select(-'...1')

# community stability
ComStab <- read_csv('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/CommunityStabilities.csv')%>%  select(-'...1') %>%
  rename(Deltabm.Tot = AUC.sum.delatbm.tot,
         Recovery=Recov,
         Resistance=Resist)

#### Response Diversity ####

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

#### Community Stability and RD ####

##### Plot - Sum Community Stability ####
AllStab <- merge(realised.pert,ComStab, by = c('caseID', 'resp.cat')) %>%
  filter(!str_detect(RD.metric,'abs' ))
str(AllStab)

#label grid
labeller <- c(mean_spp_deltabm = 'Mean Realised Response', RD_diss = 'Response Dissimilarity', RD_div = 'Response Divergence')

P_Fig4a <- AllStab %>%
  filter(!str_detect(RD.metric,'abs' ) & RD.metric == 'mean_spp_deltabm' ) %>%
  ggplot(aes(x = RD.value, y = Deltabm.Tot))+
  geom_hline(yintercept = 0)+  #  stat_poly_eq()+
  geom_point(alpha = 0.8, size = 2, color = '#F8766D')+
  labs(x = 'Mean Realised Response', y = 'OEV')+
 # facet_wrap(~RD.metric, scale = 'free_x', labeller = labeller(RD.metric = labeller))+
  theme_bw()+
  theme(axis.title.y=element_text(size=12, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=10,face="plain",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=12,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=10,face="plain",colour="black"))+
  theme(legend.position = 'none')

P_Fig4a


P_Fig4b <- AllStab %>%
  filter(!str_detect(RD.metric,'abs' ) & RD.metric == 'RD_diss') %>%
  ggplot(aes(x = RD.value, y = Deltabm.Tot))+
  geom_hline(yintercept = 0)+  #  stat_poly_eq()+
  geom_point(alpha = 0.8, size = 2, color = '#00BA38')+
  labs(x = 'Realised Response Dissimilarity', y = 'OEV')+
 # facet_wrap(~RD.metric, scale = 'free_x', labeller = labeller(RD.metric = labeller))+
  theme_bw()+
  theme(axis.title.y=element_text(size=12, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=10,face="plain",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=12,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=10,face="plain",colour="black"))+
  theme(legend.position = 'none')

P_Fig4b

P_Fig4c <- AllStab %>%
  filter(!str_detect(RD.metric,'abs' ) & RD.metric == 'RD_div') %>%
  ggplot(aes(x = RD.value, y = Deltabm.Tot))+
  geom_hline(yintercept = 0)+  #  stat_poly_eq()+
  geom_point(alpha = 0.8, size = 2, color = '#619CFF')+
  labs(x = 'Realised Response Divergence', y = 'OEV')+
  #facet_wrap(~RD.metric, scale = 'free_x', labeller = labeller(RD.metric = labeller))+
  theme_bw()+
  theme(axis.title.y=element_text(size=12, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=10,face="plain",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=12,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=10,face="plain",colour="black"))+
  theme(legend.position = 'none')
P_Fig4c

plot_grid(P_Fig4a, P_Fig4b, P_Fig4c, ncol = 3, labels = c('(a)', '(b)', '(c)'))
ggsave(plot = last_plot(), file = here('MetaMultistab/output/Fig4_RealisedResponses_OEV.png'), width = 10, height = 3.5)


#### START Meta-Analysis ####

#libraries we need
library(metafor)
library(broom)

## look at Effect sizes 
summary(AllStab)
names(AllStab)

metadata <- AllStab %>%
  filter(!str_detect(RD.metric, 'abs'))%>%
  spread(key = RD.metric, value = RD.value) 

setdiff(AllStab$caseID,metadata$caseID)
hist(metadata$OEV)
unique(metadata$studyID)

metadata <- filter(metadata, resp.cat !=  "contribution to production")

# unweighted MA requires column 1
metadata$unweighted<-1
names(metadata)



#m0
m0<-rma.mv(Deltabm.Tot,unweighted,
                      mods = ~mean_spp_deltabm+RD_diss+RD_div,
                      random = ~ 1 | caseID,
                      method="REML",data=metadata)
summary(m0) 

ModelResults <- tibble(estimate = m0$b, ci.lb = m0$ci.lb, ci.ub = m0$ci.ub, pvalue = as.numeric(m0$pval), mods = c('intercept','Mean response', 'RD dissimilarity', 'RD divergence'))
str(ModelResults)
ModelResults$mods[ModelResults$mods == 'RD divergence'] <- 'Divergence'
ModelResults$mods[ModelResults$mods == 'RD dissimilarity'] <- 'Dissimilarity'
ModelResults$mods[ModelResults$mods == 'Mean response'] <- 'Mean Response'

ModelResults$mods<- factor(ModelResults$mods, levels = c( 'Divergence','Dissimilarity','Mean Response','intercept' ))


ModelResults %>%
  mutate(p.value = as.numeric(paste(ifelse(pvalue <0.05, 1, 0.95)) ))%>%
ggplot(., aes(x = estimate, y = mods, color = mods))+
  geom_vline(xintercept = 0)+
  geom_point(size = 4 )+
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), width = .2)+
  annotate("text", y = "Mean Response", x = 1.1, label = "*", size = 8, color = "#22292F") + 
  scale_color_manual(values = c('#619CFF','#00BA38','#F8766D','#3b3b3b'))+
  labs(y = '')+
  theme_bw()+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,colour="black",angle=0,hjust=0.4),
        axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,colour="black"),
        panel.border=element_rect(colour="black",linewidth=1.5),
        legend.position = 'none')
ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/Metamultistab/output/Forestplot_gg.pdf'), width = 6, height = 4)




#### Appendix: Stability Metrics ####

Metrics <- AllStab %>%
  filter(!str_detect(RD.metric,'abs' ) & RD.metric != 'var_spp_deltabm')

Metrics$RD.metric<- factor(Metrics$RD.metric, levels = c( 'mean_spp_deltabm','RD_diss','RD_div'))

labeller <- c(mean_spp_deltabm = 'Mean Realised Response', RD_diss = 'Realised Response Dissimilarity', RD_div = 'Realised Response Divergence')

p2<-Metrics %>%
  ggplot(., aes ( y = Resistance, x = RD.value, color = RD.metric))+
  scale_color_manual(values = c('#F8766D','#00BA38','#619CFF'))+
  labs(x='')+
  geom_hline(yintercept = 0)+
  geom_point(size = 2, alpha = 0.8)+
  facet_wrap(~RD.metric, scales='free_x', labeller = labeller(RD.metric = labeller))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = 'none')
p2
#ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/ResponseTraits_Resist_sum.png'), width =  6, height = 5)

p3<-Metrics %>%
  ggplot(., aes ( y = Recovery, x = RD.value, color = RD.metric))+
  scale_color_manual(values = c('#F8766D','#00BA38','#619CFF'))+
  labs(x='')+
  geom_hline(yintercept = 0)+
  geom_point(size = 2, alpha = 0.8)+
  facet_wrap(~RD.metric, scales='free_x', labeller = labeller(RD.metric = labeller))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = 'none')
p3


p4<-Metrics %>%
  ggplot(., aes ( y = CV, x = RD.value, color = RD.metric))+
  scale_color_manual(values = c('#F8766D','#00BA38','#619CFF'))+
  labs(x='')+
  geom_hline(yintercept = 0)+
  geom_point(size = 2, alpha = 0.8)+
  facet_wrap(~RD.metric, scales='free_x', labeller = labeller(RD.metric = labeller))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = 'none')
p4

plot_grid(p2,p3,p4, ncol = 1, labels = c('(a)', '(b)', '(c)'))
ggsave(plot = last_plot(), file= here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/Appendix_FigS_SumStabilityMetric.png'), width = 7, height = 8)



