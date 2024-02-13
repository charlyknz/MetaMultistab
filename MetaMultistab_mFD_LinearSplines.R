#### R script for meta analysis ####

rm(list=ls())
graphics.off()

library(tidyverse)
library(readxl)
library(MESS)
library(here)
library(cowplot)
library(GGally)
library(ggpubr)

#### import data ####
setwd("~/Desktop/phD/Meta_Multistab")
study <- read_excel("~/Desktop/phD/Meta_Multistab/MetaMultistab/Multistab_species_data_mfd.xlsx") %>%
  select(-c(35:51)) %>%
  rename(caseID = paper)

rawData <- read_excel("~/Desktop/phD/Meta_Multistab/MetaMultistab/Multistab_species_data_mfd.xlsx", 
                      sheet = "species data")

allData <- study%>%
  select(-func, -resp.cat,-resp) %>%
  merge(., rawData,by = c('caseID'))%>%
  filter(spec.inf %in% c('species', 'taxa')) %>%
  drop_na(caseID)

which(is.na(allData$Con.M))

communityStab <- read_excel("~/Desktop/phD/Meta_Multistab/MetaMultistab/Multistab_species_data_mfd.xlsx", 
                            sheet = "communityFunction") %>%
  select(caseID, totRR, resp, DAY, RD)

#### test if all data are imported correctly ####
unique(rawData$caseID)
unique(allData$studyID)
unique(allData$caseID)


setdiff( study$caseID,allData$caseID) 
# CK028 sum of difference
# CK035 press disturbance and 1 timepoint 

# look at data
unique(allData$species_specification)
allData$Con.M[(allData$Con.M <0)]<-0
allData$Dist.M[(allData$Dist.M <0)]<-0

response <- allData %>%
  select(caseID, studyID, spec.inf, comment.x,  system, lat, long, organism, duration, differentiation,dist.cat,open, species, species_specification,func,resp, resp.cat,DAY, RD,Con.M, Dist.M,Con.N, Dist.N)%>%
  mutate(dummyRR = Con.M + Dist.M) %>% 
  filter(dummyRR != 0) %>%##take out those rows where biomass is 0 in both treatment (Biomass) + control (con.bio)
  group_by(caseID, RD)%>%
  mutate(con.tot = sum(Con.M, na.rm = T),
         dist.tot = sum(Dist.M, na.rm = T))%>%
  ungroup() %>%
  mutate(dist.pi = Dist.M/dist.tot,
         con.pi = Con.M/con.tot,
         delta.pi = dist.pi - con.pi, # calculate effect sizes
         RR = (Dist.M-Con.M)/(Dist.M+Con.M),
         LRR = log(dist.tot/ con.tot),
         deltabm.tot = (dist.tot - con.tot)/(dist.tot+con.tot)) %>%
  mutate(USI = paste(caseID, species,  sep = "_"))  %>%
  filter(resp.cat != "contribution to production") %>%
  distinct(caseID, studyID, spec.inf, comment.x,  system, lat, long, organism, duration, differentiation,dist.cat,open, species, species_specification,func,resp, resp.cat,DAY, RD,Con.M, Dist.M,Con.N, Dist.N,
           deltabm.tot,LRR,RR,delta.pi,con.pi,dist.pi,USI)

response$delta.pi[is.infinite(response$delta.pi)]<-NA
response$RR[is.infinite(response$RR)]<-NA

hist(response$RR)


## check if meta analysis mods are not NA
which(is.na(response$duration))
which(is.na(response$open))
which(is.na(response$resp.cat))
which(is.na(response$dist.cat))

#### AUC Loop ####
USI <- response$USI
response <- response[order(response$RD),]
source("/Users/charlottekunze/Desktop/phD/Meta_Multistab/myauc.R")

names(response)

stab.auc <- tibble()

for(i in 1:length(USI)){
  temp<-response[response$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>3){#does the next step only if at least 3 data points are present
    AUC.RR<-myauc(temp$RD, temp$RR,  from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                absolutearea = FALSE)
    AUC.pi<-myauc(temp$RD, temp$delta.pi, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                absolutearea = FALSE)
    AUC.totRR<-myauc(temp$RD, temp$deltabm.tot, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                   absolutearea = TRUE)
    AUC.RR.spline<-myauc(temp$RD, temp$RR,  from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                       type = c("spline"),absolutearea = FALSE)
    AUC.pi.spline<-myauc(temp$RD, temp$delta.pi, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                       type = c("spline"),absolutearea = FALSE)
    AUC.totRR.spline<-myauc(temp$RD, temp$deltabm.tot, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                          type = c("spline"),absolutearea = TRUE)
    con.pi <- mean(temp$con.pi, na.rm = T)
    mean.delta.pi <- mean(temp$delta.pi, na.rm = T)
    mean.RR <- mean(temp$RR, na.rm = T)
    mean.totRR <- mean(temp$totRR, na.rm = T)
    sum.con <- sum(temp$Con.M)
    stab.auc<-rbind(stab.auc,
                    tibble(temp[1,c(1:17)],
                           AUC.RR = AUC.RR[[1]],
                           AUC.pi = AUC.pi[[1]],
                           AUC.RR.spline = AUC.RR.spline[[1]],
                           AUC.pi.spline = AUC.pi.spline[[1]],
                           con.pi,
                           sum.con,
                           mean.RR, mean.delta.pi,mean.totRR,
                           AUC.totRR = AUC.totRR[[1]],
                           AUC.totRR.spline= AUC.totRR.spline[[1]]))
    rm(temp)
  }
}


#### remove duplicates ####
unique(stab.auc$caseID)

str(stab.auc)

data.plot <- stab.auc %>%
  distinct(caseID, studyID,system, organism, duration, dist.cat, open, species, resp.cat, sum.con,AUC.RR, AUC.RR.spline, AUC.pi.spline, AUC.pi, AUC.totRR, con.pi,mean.totRR, mean.delta.pi,mean.RR)%>%
  filter(resp.cat %in% c('abundance', 'biomass') ) %>%
  ungroup() %>%
  drop_na(AUC.pi)
names(data.plot)
hist(data.plot$AUC.RR)
hist(data.plot$AUC.RR.spline)


#### spline and linear methods ###
linearReg <- lm(data.plot$AUC.RR.spline~data.plot$AUC.RR)
summary(linearReg)

linearRegPi <- lm(data.plot$AUC.pi.spline~data.plot$AUC.pi)
summary(linearRegPi)

ggplot(data.plot, aes(y = AUC.pi, x = AUC.pi.spline))+
  geom_point(size = 1.5)+
  labs(x= 'AUC.pi natural splines', y = 'AUC.pi linear splines')+
  scale_x_continuous(limits = c(-3,4))+
  scale_y_continuous(limits = c(-3,4))+
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()
#ggsave(plot = last_plot(), file = 'CorrelationAUC.pi_LinearNatural.png')

ggplot(data.plot, aes(y = AUC.RR, x = AUC.RR.spline))+
  geom_point(size = 1.5)+
  labs(x = 'AUC.RR natural splines', y= 'AUC.RR linear splines')+
  scale_x_continuous(limits = c(-10,5))+
  scale_y_continuous(limits = c(-10,5))+
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()
#ggsave(plot = last_plot(), file = 'CorrelationAUC.RR_LinearNatural.png')



#1
raw <-ggplot(data.plot, aes(x=AUC.pi,y=AUC.RR, col=system))+
  geom_hline(yintercept=0, colour="grey")+
  geom_vline(xintercept=0, colour="grey")+
  geom_point(alpha=0.4, size = 2)+
  scale_shape_manual(values=c(16, 17,1,2,5))+
  theme_bw()+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(colour="black",linewidth =1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",linewidth=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))#+
# geom_smooth(col="black")
raw
ggsave(plot = raw, file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/AUC.RRdeltaPi.png'), width = 4, height = 4)

#1
raw1 <-ggplot(data.plot, aes(x=mean.delta.pi,y=mean.RR, col=system))+
  geom_hline(yintercept=0, colour="grey")+
  geom_vline(xintercept=0, colour="grey")+
  geom_point(alpha=0.4, size = 2)+
  scale_shape_manual(values=c(16, 17,1,2,5))+
  theme_bw()+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(colour="black",linewidth =1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",linewidth=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))#+
# geom_smooth(col="black")
raw1
#ggsave(plot = raw1, file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/MeanRRdeltaPi.png'))



#### COUnt SECTORS ####
data.plot$Sector<- NA
data.plot$Sector[data.plot$mean.delta.pi>0&data.plot$mean.RR>0]<-1
data.plot$Sector[data.plot$mean.delta.pi>0&data.plot$mean.RR<0]<-2
data.plot$Sector[data.plot$mean.delta.pi<0&data.plot$mean.RR<0]<-3
data.plot$Sector[data.plot$mean.delta.pi<0&data.plot$mean.RR>0]<-4
data.plot$Sector[data.plot$mean.delta.pi== 0 & data.plot$mean.RR==0] <-0


sector.count<- data.plot%>%
  group_by( Sector) %>%
  count()%>%
  rename(N = n) %>%
  drop_na(Sector) %>%
  ungroup()%>%
  mutate(sum = sum(N, na.rm = T),
         relN = (N/sum)*100)
sector.count


ggplot(data.plot,aes(AUC.pi, AUC.RR, color=resp.cat )) +
  geom_point(alpha = 0.4,  size = 3) +
  geom_vline(xintercept = 0, alpha = 0.5) +                                      
  geom_hline(yintercept = 0, alpha = 0.5) +
  facet_wrap(~system, scales = 'free') +
  labs(x = 'Relative Contribution to Stability', y = "Absolute Contribution to Stability", color = 'System', shape = 'System') +  
  theme_bw()+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(colour="black",linewidth=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))
ggsave(plot = last_plot(), file = here('~/Desktop/PhD/Meta_Multistab/MetaMultistab/output/AUCpiAUCRR_system.png'), width = 8, height = 5)


hist(data.plot$AUC.RR)

#### dominance ####
ggplot(data.plot,aes(con.pi, AUC.RR, shape = system, color=resp.cat )) +
  geom_point(alpha = 0.4,  size = 3) +
  geom_vline(xintercept = 0, alpha = 0.5) +                                      
  geom_hline(yintercept = 0, alpha = 0.5) +
  #facet_wrap(~system, scales = 'free') +
  labs(x = 'relative dominance', y = "Absolute Contribution to Stability", color = 'System', shape = 'System') +  
  theme_bw()+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(colour="black",linewidth=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))
ggsave(plot = last_plot(), file = here('~/Desktop/PhD/Meta_Multistab/MetaMultistab/output/CONpiAUCRR_system.png'), width = 8, height = 5)


#### dominance ####
ggplot(data.plot,aes(con.pi, AUC.pi, shape = system, color=resp.cat )) +
  geom_point(alpha = 0.4,  size = 3) +
  geom_vline(xintercept = 0, alpha = 0.5) +                                      
  geom_hline(yintercept = 0, alpha = 0.5) +
  #facet_wrap(~system, scales = 'free') +
  labs(x = 'relative dominance', y = "Relative Contribution to Stability", color = 'System', shape = 'System') +  
  theme_bw()+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(colour="black",linewidth=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))
ggsave(plot = last_plot(), file = here('~/Desktop/PhD/Meta_Multistab/MetaMultistab/output/CONpiAUCpi_system.png'), width = 8, height = 5)

hist(data.plot$AUC.RR)



hist(data.plot$AUC.RR)

#### richness ####

SR <- data.plot %>%
  count(caseID)

data.plot %>%
  left_join(., SR) %>%
  gather(AUC.RR, AUC.pi, key = 'metric', value = 'value') %>%
ggplot(.,aes(n, value, shape = system, color=resp.cat )) +
  geom_point(alpha = 0.4,  size = 3) +
  geom_vline(xintercept = 0, alpha = 0.5) +                                      
  geom_hline(yintercept = 0, alpha = 0.5) +
  facet_wrap(~metric, scales = 'free') +
  labs(x = 'S', color = 'System', shape = 'System') +  
  theme_bw()+
 # scale_x_continuous(trans='log')+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(colour="black",linewidth=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))
ggsave(plot = last_plot(), file = here('~/Desktop/PhD/Meta_Multistab/MetaMultistab/output/SRAUCs.png'), width = 8, height = 5)

SR <- data.plot %>%
  count(caseID)

data.plot %>%
  left_join(., SR) %>%
  gather(AUC.RR, AUC.pi, key = 'metric', value = 'value') %>%
  ggplot(.,aes(con.pi, value, shape = system, color=n )) +
  geom_point(alpha = 0.4,  size = 3) +
  geom_vline(xintercept = 0, alpha = 0.5) +                                      
  geom_hline(yintercept = 0, alpha = 0.5) +
  facet_wrap(~metric, scales = 'free') +
  labs(x = 'relative dominance', S = 'System', shape = 'System') +  
  theme_bw()+
  # scale_x_continuous(trans='log')+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(colour="black",linewidth=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))
ggsave(plot = last_plot(), file = here('~/Desktop/PhD/Meta_Multistab/MetaMultistab/output/SR_conPI_AUCs.png'), width = 8, height = 5)

#### Community Stability from Meta-Analysis ####
unique(communityStab$resp)

communityStab1<- communityStab %>%
  filter(!resp %in% c('production', 'respiration')) %>%
  drop_na(totRR) %>%
  rename(resp.cat = resp)


communityStab1$resp.cat[communityStab1$resp.cat == 'cover'] <-'biomass'

## create USI to run loop
USIc <- communityStab1$caseID
names(communityStab1)

stab.auc.c <- data.frame()
for(i in 1:length(USIc)){
  temp<-communityStab1[communityStab1$caseID==USIc[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>3){#does the next step only if at least 3 data points are present
    AUC.tot.RR<-auc(temp$RD, temp$totRR, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                    type = c("linear"),absolutearea = TRUE)
    AUC.totRR.spline <- auc(temp$RD, temp$totRR, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                            type = c("spline"),absolutearea = TRUE)
    mean.tot.RR <- mean(temp$totRR, na.rm = T)
    nrow <- nrow(temp)
    stab.auc.c<-rbind(stab.auc.c,
                      data.frame(temp[1,c(1,3)],
                                 AUC.tot.RR,mean.tot.RR,
                                 AUC.totRR.spline, nrow))
    rm(temp)
  }
}

summary(stab.auc.c)

communityAUC <- distinct(stab.auc.c, caseID, resp.cat, AUC.tot.RR, mean.tot.RR,nrow, AUC.totRR.spline)
ggplot(communityAUC, aes(x = AUC.tot.RR, y = AUC.totRR.spline))+
  geom_point(size = 1.5)+
  labs(x = 'AUC natural splines', 'AUC linear splines')+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,5))+
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()
#ggsave(plot = last_plot(), file = 'CorrelationAUC_OEVLinearNatural.png')

lm1 <-lm(communityAUC$AUC.totRR.spline~communityAUC$AUC.tot.RR)
summary(lm1)



## join Meta-Analysis stability values with sum of individual species if communitystab is not present ##
all.stab.auc <- left_join(data.plot,communityAUC,  by = c('caseID','resp.cat')) %>%
  mutate(OEV = ifelse(is.na(AUC.tot.RR), AUC.totRR, AUC.tot.RR ),
         meanTotalRR = ifelse(is.na(mean.tot.RR), mean.totRR, mean.tot.RR ),) %>%
  distinct(caseID,studyID, organism, system,duration, dist.cat, open, resp.cat, OEV,meanTotalRR)
which(is.na(all.stab.auc$OEV))

#### Linear and natural splines difference ####
str(communityAUC)
comparison.df <- communityAUC%>%
  merge(., data.plot, by = c('caseID','resp.cat')) %>%
  mutate(diff = AUC.totRR.spline - AUC.tot.RR ) 

which(is.na(comparison.df$diff))
n <- count(comparison.df, nrow)
unique(comparison.df$caseID)

ggplot(comparison.df, aes(x = nrow, y = diff))+
  geom_hline(yintercept = 1, color = 'grey', linetype = 'dashed', linewidth = 0.5)+
  geom_hline(yintercept = 0.5, color = 'grey', linetype = 'dashed', linewidth = 0.5)+
  geom_point()+
  labs(x = 'observations', y = 'Difference OEV natural - linear splines')+
  theme_bw()
#ggsave(plot=last_plot(), file = 'DifferenceLinearNatural.png')

#### ResponseDiversity ####
names(data.plot)
source(here("~/Desktop/phD/response-diversity-pulse-pert/R/Ross_et_al_functions.R"))

igr.pert <- data.plot %>%
  group_by(caseID) %>%
   mutate(
         mean_igr_effect = mean(AUC.RR),
         var_igr_effect = var(AUC.RR),
         RD_diss = resp_div(AUC.RR, sign_sens = FALSE),
         RD_div = resp_div(AUC.RR, sign_sens = TRUE)) %>%
  ungroup()%>%
  gather(mean_igr_effect,var_igr_effect,RD_diss,RD_div, key = 'RD.metric', value = 'RD.value') %>%
  left_join(., all.stab.auc) %>%
  distinct(RD.metric, RD.value, OEV, meanTotalRR,caseID, studyID, system, organism, duration, dist.cat, open,
           resp.cat)

ggplot(igr.pert, aes ( x = RD.value, y = OEV))+
  geom_point()+
  facet_wrap(~RD.metric, scales = 'free_x')+
  theme_bw()
ggsave(plot=last_plot(), file = here('MetaMultistab/output/MeanIGReffect.png'), width =  6, height = 4)





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
  select(population, mean.delta.pi, mean.RR)%>%
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
  faxes           = c('mean.RR','mean.delta.pi'),
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
  left_join(., all.stab.auc, by = c('caseID')) 
which(is.na(Fdiv$OEV))

Fdiv$organism[Fdiv$organism == 'Cladoceran'] <- 'zooplankton'
Fdiv$organism[Fdiv$organism == 'invertebrate'] <- 'invertebrates'
Fdiv$organism[Fdiv$organism == 'macroinvertebrate'] <- 'invertebrates'

unique(Fdiv$organism)
Fdiv$OrganismType <- Fdiv$organism
Fdiv$OrganismType[Fdiv$OrganismType == 'macroalgae'] <- 'primaryproducer'
Fdiv$OrganismType[Fdiv$OrganismType == 'phytoplankton'] <- 'primaryproducer'

Fdiv$OrganismType[Fdiv$OrganismType == 'periphyton'] <- 'primaryproducer'
Fdiv$OrganismType[Fdiv$OrganismType == 'macrophytes'] <- 'primaryproducer'
Fdiv$OrganismType[Fdiv$OrganismType == 'plant'] <- 'primaryproducer'
Fdiv$OrganismType[Fdiv$OrganismType == 'zooplankton'] <- 'herbivore'
Fdiv$OrganismType[Fdiv$OrganismType == 'vertebrates'] <- 'predator'
Fdiv$OrganismType[Fdiv$OrganismType == 'fish'] <- 'predator'
Fdiv$OrganismType[Fdiv$OrganismType == 'invertebrates'] <- 'predator'
Fdiv$OrganismType[Fdiv$OrganismType == 'invertebrate'] <- 'predator'
Fdiv$OrganismType[Fdiv$OrganismType == 'macroinvertebrate'] <- 'predator'


str(Fdiv)
write.csv(Fdiv, file = 'Output/StabFdiv_n4.csv')

### Plot FDiv Indices ~ Stability ###

Fdiv %>%
  select(caseID, fdis, fric, feve, sp_richn,studyID, organism, dist.cat, OEV, resp.cat,system,OrganismType) %>%
  gather(c(fdis, fric, feve), key = 'Indices', value = 'IndexValue') %>%
  ggplot(., aes(x = IndexValue, y = OEV))+
  geom_point(size = 2, alpha = .4)+
  labs(x = 'F Dispersion', y = 'OEV')+
  facet_grid(~Indices, scales = 'free_x') +
  theme_bw()+
  theme(text = element_text(size=rel(4)),
        strip.text.x = element_text(size=rel(4)))
ggsave(last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/Fdiv_overview.png'), width = 8, height = 5)
Fdiv %>%
  select(caseID, fdis, fric, feve, sp_richn,studyID, organism, dist.cat, OEV, resp.cat,system,OrganismType) %>%
  gather(c(fdis, fric, feve), key = 'Indices', value = 'IndexValue') %>%
  ggplot(., aes(x = IndexValue, y = OEV))+
  geom_hline(yintercept = 0)+
  geom_smooth(method='lm')+
  geom_point(size = 2, alpha = .4)+
  labs(x = 'F Dispersion', y = 'OEV')+
  facet_grid(~Indices, scales = 'free_x') +
  theme_bw()+
  theme(legend.position = 'bottom')
ggsave(last_plot(), file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/Fdiv.png'), width = 10, height = 5)

fdis<-ggplot(Fdiv, aes(x = fdis, y = OEV))+
  geom_hline(yintercept = 0)+
  geom_point(size = 2, alpha = .4)+
  labs(x = 'F Dispersion', y = 'OEV')+
  theme_bw()+
  geom_smooth(method='lm')+
  theme(legend.position = 'bottom')

fric<-ggplot(Fdiv, aes(x = fric, y = OEV))+
  geom_hline(yintercept = 0)+
  geom_point(size = 2, alpha = .4)+
  geom_smooth(method='lm')+
  labs(x = 'F Richness', y = 'OEV')+
  theme_bw()+
  theme(legend.position = 'bottom')

feve<-ggplot(Fdiv, aes(x = feve, y = OEV))+
  geom_hline(yintercept = 0)+
  geom_point(size = 2, alpha = .4)+
  geom_smooth(method='lm')+
  labs(x = 'F Evenness', y = 'OEV')+
  theme_bw()+
  theme(legend.position = 'bottom')

sp_richn<-ggplot(Fdiv, aes(x = sp_richn, y = OEV))+
  geom_hline(yintercept = 0)+
  geom_point(size = 2, alpha = .4)+
  geom_smooth(method='lm')+
  labs(x = 'S Richness', y = 'OEV')+
  theme_bw()+
  theme(legend.position = 'bottom')

allPlot <- plot_grid(fric+ theme(legend.position="none"),feve+ theme(legend.position="none"),fdis+ theme(legend.position="none"), sp_richn+ theme(legend.position="none"))
legend_b <- get_legend(fric + theme(legend.position="right"))

plot_grid(allPlot, legend_b,rel_widths = c(1,0.3))

ggsave(last_plot(), file = here('~/Desktop/phD/Meta_Multistab/FDiv_LinFdivAdj.png'), width = 14, height = 10)


plot1 <- ggscatter(Fdiv, x = 'fric', y='OEV', cor.coef = T, add= 'reg.line')
plot2 <- ggscatter(Fdiv, x = 'feve', y='OEV', cor.coef = T, add= 'reg.line')
plot3 <- ggscatter(Fdiv, x = 'fdis', y='OEV', cor.coef = T, add= 'reg.line')

plot1+plot2+plot3

#### START Meta-Analysis ####

# import data 
Fdiv <- read.csv('Output/StabFdiv.csv')

#libraries we need
library(metafor)
library(psych)
library(grid)
library(gridExtra)
library(reshape2)
library(agricolae)
library(ggExtra)
library(metafor)
library(cowplot)
library(broom)
library(GGally)


## look at Effect sizes 

summary(Fdiv)
names(Fdiv)

metadata <- distinct(Fdiv, caseID ,studyID,duration, open, organism, OrganismType, system,dist.cat,resp.cat,OEV, fdis, fric, feve, sp_richn)
hist(metadata$OEV)
unique(metadata$caseID)

metadata <- filter(metadata, resp.cat !=  "contribution to production")


# unweighted MA requires column 1
metadata$unweighted<-1

test.complete<-rma.mv(OEV,unweighted,
                      mods = ~fdis+fric+feve+resp.cat+dist.cat+system+open+organism+OrganismType+duration,
                      random = ~ 1 | caseID,
                      method="REML",data=metadata)
summary(test.complete)
complete.mod<-tidy(summary(test.complete))
complete.mod$k<-test.complete$k
complete.mod$AIC<-test.complete$fit.stats$REML[3]

