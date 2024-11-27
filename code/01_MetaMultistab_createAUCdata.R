#### R script for meta analysis ####

setwd("~/Desktop/phD/Meta_Multistab")

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
  select(caseID, totRR, resp, DAY, RD, RESIST, RESIL, RECOV)

# test if all data are imported correctly
unique(rawData$caseID)
unique(allData$studyID)
unique(allData$caseID)


setdiff( study$caseID,allData$caseID) 
# CK028 sum of difference
# CK035 press disturbance and 1 timepoint 

# look at data
unique(allData$species_specification)

#remove negative biomass and abundance entries!
allData$Con.M[(allData$Con.M <0)]<-0
allData$Dist.M[(allData$Dist.M <0)]<-0

response <- allData %>%
  select(caseID, studyID, spec.inf, comment.x,  system, lat, long, organism, duration, differentiation,dist.cat,open, 
         species, species_specification,func,resp, resp.cat,DAY, RD,Con.M, Dist.M,Con.N, Dist.N, Dist.SD, Con.SD)%>%
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
           deltabm.tot,LRR,RR,delta.pi,con.pi,dist.tot, con.tot,dist.pi,USI)


#remove infinite values for species relative proportion delta.pi and absolute biomass response ratio rr
response$delta.pi[is.infinite(response$delta.pi)]<-NA
response$RR[is.infinite(response$RR)]<-NA
response$RR[response$RR == 'NaN']<-NA
hist(response$RR)


## check if meta analysis mods are not NA
which(is.na(response$duration))
which(is.na(response$open))
which(is.na(response$resp.cat))
which(is.na(response$dist.cat))


#### AUC Loop species contributions ####
USI <- response$USI #unique identifier

#order after time steps
response <- response[order(response$RD),]

#auc function which allows us to plot predictions
source("/Users/charlottekunze/Desktop/phD/Meta_Multistab/myauc.R")

names(response)

#create empty df
stab.auc <- tibble()

for(i in 1:length(USI)){
  temp<-response[response$USI==USI[i], ]#creates a temporary data frame for each case
   if(dim(temp)[1]>3){#does the next step only if at least 3 data points are present
    AUC.RR<-auc(temp$RD, temp$RR,  from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                absolutearea = FALSE)
    AUC.pi<-auc(temp$RD, temp$delta.pi, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                absolutearea = FALSE)
    con.pi <- mean(temp$con.pi, na.rm = T)
    mean.delta.pi <- mean(temp$delta.pi, na.rm = T)
    mean.RR <- mean(temp$RR, na.rm = T)
    sum.con <- sum(temp$Con.M)
     stab.auc<-rbind(stab.auc,
                    tibble(temp[1,c(1:17)],
                           AUC.RR ,
                           AUC.pi ,
                           con.pi,
                           sum.con,
                           mean.RR, 
                           mean.delta.pi
                          ))
    rm(temp)
  }
}

### remove duplicates ###
unique(stab.auc$caseID)
str(stab.auc)

data.plot <- stab.auc %>%
  distinct(caseID, studyID,system, lat, organism, duration, dist.cat, open, species, resp.cat, sum.con,AUC.RR,  AUC.pi, con.pi, mean.delta.pi,mean.RR)%>%
  filter(resp.cat %in% c('abundance', 'biomass') ) %>%
  ungroup() %>%
  drop_na(AUC.pi)%>% 
  #remove studies with less than 2 species # Note: sometimes species get removed during the AUC loop
  group_by(caseID) %>% 
  filter(n() >1)
#rm(stab.auc)

names(data.plot)
hist(data.plot$AUC.RR)
hist(data.plot$AUC.pi)

#write.csv(data.plot, file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/SpeciesStabilities.csv'))


### For resistance we have to slice the first entry after the start community ##
# as sometimes we dont have information on the start community, we will split the two df 
resist <- com.response%>%
  group_by(caseID) %>%
  filter('Start' %in% Stab.metric) %>%
  arrange(caseID, RD) %>%
  group_by(caseID) %>%
  slice(2)

resistance <- com.response%>%
  group_by(caseID) %>%
  filter(!'Start' %in% Stab.metric) %>%
  arrange(caseID, RD) %>%
  group_by(caseID) %>%
  slice(1) %>%
  bind_rows(., resist) %>%
  rename(Resist_LRR = LRR, 
         Resist = deltabm.tot) %>%
  select(caseID, resp.cat, Resist, Resist_LRR)
resistance$Resist_LRR[is.infinite(resistance$Resist_LRR)]<-NA

#create USI
caseID <- com.response$caseID
unique(caseID)

com.stab <- tibble()

for(i in 1:length(caseID)){
  temp<-com.response[com.response$caseID==caseID[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>3){#does the next step only if at least 3 data points are present
    OEV<-auc(temp$RD, temp$deltabm.tot, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE), 
                 absolutearea = TRUE)
    AUC.sum.delatbm.tot <-auc(temp$RD, temp$deltabm.tot, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                 type = c("linear"),absolutearea = FALSE)
    CV<- mean(temp$deltabm.tot, na.rm = T)/sd(temp$deltabm.tot, na.rm = T) # coefficient of variation
    com.stab<-rbind(com.stab,
                    tibble(temp[1,c(1:12)],
                           OEV ,CV,
                           AUC.sum.delatbm.tot))
    rm(temp)
  }
}

str(com.stab)

#recovery
recov <- com.response %>%
  filter(RD == 1) %>%
  rename(Recov = deltabm.tot) %>%
  distinct(caseID, resp.cat,Recov)

com.stab.sum <- com.stab %>%
  distinct(caseID, resp.cat, AUC.sum.delatbm.tot,OEV,CV) %>%
  left_join(., recov) %>%
  left_join(., resistance)
names(com.stab)

#rm(com.stab)

#write.csv(com.stab.sum, file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/CommunityStabilities.csv'))

### variance and community stability ###

spp.var <- allData%>%
  distinct(caseID, studyID,  species, resp.cat, RD,Con.M, Dist.M,Con.N, Dist.N, Dist.SD, Con.SD)%>%
  mutate(dummyRR = Con.M + Dist.M) %>% 
  filter(dummyRR != 0) %>%##take out those rows where biomass is 0 in both treatment (Biomass) + control (con.bio)
  group_by(caseID, resp.cat, species)%>%
  reframe(Mean.Dist = mean(Dist.M,na.rm = T),
            sd.Dist = sd(Dist.M, na.rm = T),
            Mean.Con = mean(Con.M,na.rm = T),
            sd.Con = sd(Con.M, na.rm = T),
         CV.Dist= (Mean.Dist/as.numeric(sd.Dist)),
         CV.con= (Mean.Con/as.numeric(sd.Con)),
         Dist.N = Dist.N, 
         Con.N = Con.N)%>%
  mutate(lnCVR= (log(CV.Dist/CV.con) + (1/(2*(Dist.N-1)))-(1/(2*(Dist.N-1)))),
         lnCV = NA) %>% #lnCV =  ln sd  ln mean 1 +(1/(2n-1))
  drop_na(lnCVR)%>%
  filter(!is.infinite(lnCVR))
  
varStab <- com.stab %>%
  distinct(caseID, resp.cat, system, organism,AUC.sum.delatbm.tot,OEV,CV) %>%
  left_join(., spp.var)
  
varStab$organisms <- varStab$organism
varStab$organism[varStab$organism %in% c('macroinvertebrate', 'invertebrates','fish','vertebrates')] <- 'predator'
varStab$organism[varStab$organism%in%c('plant','macroalgae','macrophytes','periphyton')] <- 'primaryproducer'
varStab$organism[varStab$organism %in% c('Copepod','Cladoceran','beetle','zooplankton','meiofauna')] <- 'herbivore'

varStab %>%
   ggplot(., aes(x = lnCVR, y = AUC.sum.delatbm.tot, color = organism))+
  geom_hline(yintercept = 0)+
  geom_point(size = 2)+
  facet_grid(~organism )+
  labs(y = 'Instability (OEV)', x = 'Species variance', col = 'organism')+
  theme_bw()+
  theme(legend.position='bottom')
ggsave(plot = last_plot(), file = here('lnCVR_trophy.png'), width = 8, height = 5)


varStab %>%
  group_by(caseID, AUC.sum.delatbm.tot, system, organism) %>%
  summarise(mean = mean(lnCVR)) %>%
ggplot(., aes(x = mean, y = AUC.sum.delatbm.tot))+
  geom_hline(yintercept = 0)+
  geom_point()+
  labs(y = 'Instability (OEV)', x = 'Mean species variance')+
  theme_bw()
ggsave(plot = last_plot(), file = here('MetaMultistab/output/Mean_lnCVR_plot.png'))

##### species plots #####

#AUC for systems
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
  theme(panel.border=element_rect(colour="black",linewidth=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))
ggsave(plot = last_plot(), file = here('~/Desktop/PhD/Meta_Multistab/MetaMultistab/output/AUCpiAUCRR_system.png'), width = 8, height = 5)


#AUC for organism
unique(data.plot$organism)

data.plot$organisms <- data.plot$organism
data.plot$organism[data.plot$organism %in% c('macroinvertebrate', 'invertebrates','fish','vertebrates')] <- 'predator'
data.plot$organism[data.plot$organism%in%c('plant','macroalgae','macrophytes','periphyton')] <- 'primaryproducer'
data.plot$organism[data.plot$organism %in% c('Copepod','Cladoceran','beetle','zooplankton','meiofauna')] <- 'herbivore'

ggplot(data.plot,aes(AUC.pi, AUC.RR, color=system )) +
  geom_point(alpha = 0.4,  size = 3) +
  geom_vline(xintercept = 0, alpha = 0.5) +                                      
  geom_hline(yintercept = 0, alpha = 0.5) +
  facet_wrap(~organism, scales = 'free') +
  labs(x = 'Relative Contribution to Stability', y = "Absolute Contribution to Stability", color = 'System', shape = 'System') +  
  theme_bw()+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(colour="black",linewidth=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",linewidth=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))
ggsave(plot = last_plot(), file = here('~/Desktop/PhD/Meta_Multistab/MetaMultistab/output/AUCpiAUCRR_organismtype.png'), width = 8, height = 5)


# Mean 
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



#### Count occurrences in SECTORS ####
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



##### Species Dominance #####
p1<-ggplot(data.plot,aes(con.pi, AUC.RR, shape = system, color=organism )) +
  geom_point(alpha = 0.4,  size = 3) +
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


p2<-ggplot(data.plot,aes(con.pi, AUC.pi, shape = system, color=organism )) +
  geom_point(alpha = 0.4,  size = 3) +
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

plot_grid(p1,p2,ncol = 2)
ggsave(plot=last_plot(), file = here('~/Desktop/PhD/Meta_Multistab/MetaMultistab/output/CONpiAUC_organism.png'), width = 10, height = 5)


##### Species Richness #####

SR <- data.plot %>%
  count(caseID) %>%
  left_join(., com.stab.sum)

ggscatter(SR, x = 'n', y = 'OEV', add = 'reg.line', cor.coef = T)

names(com.stab.sum)
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
#ggsave(plot = last_plot(), file = here('~/Desktop/PhD/Meta_Multistab/MetaMultistab/output/SRAUCs.png'), width = 8, height = 5)


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



#### AUC Loop Community Stability MA data ####
unique(communityStab$resp)

communityStab1<- communityStab %>%
  filter(!resp %in% c('production', 'respiration')) %>%
  drop_na(totRR) %>%
  rename(resp.cat = resp)%>%
  mutate(Stab.metric = paste(ifelse(RD == 1, 'Recovery', ifelse(RD==0, 'Start', ''))))

communityStab1$resp.cat[communityStab1$resp.cat == 'cover'] <-'biomass'

## create USI to run loop
USIc <- communityStab1$caseID
names(communityStab1)

#empty df
com.stab.MA <- data.frame()

for(i in 1:length(USIc)){
  temp<-communityStab1[communityStab1$caseID==USIc[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>3){#does the next step only if at least 3 data points are present
    OEV.MA<-auc(temp$RD, temp$totRR, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                    type = c("linear"),absolutearea = TRUE)
    AUC.delatbm.tot.MA  <-auc(temp$RD, temp$totRR, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                               type = c("linear"),absolutearea = FALSE)
    CV.MA<- mean(temp$totRR, na.rm = T)/sd(temp$totRR, na.rm = T) # coefficient of variation
    com.stab.MA<-rbind(com.stab.MA,
                      data.frame(temp[1,c(1,3)],
                                 OEV.MA,CV.MA,AUC.delatbm.tot.MA))
    rm(temp)
  }
}
str(com.stab.MA)

### For resistance we have to slice the first entry after the start community ##
# as sometimes we dont have information on the start community, we will split the two df 
resist1 <- communityStab1%>%
  group_by(caseID) %>%
  filter('Start' %in% Stab.metric) %>%
  arrange(caseID, RD) %>%
  group_by(caseID) %>%
  slice(2)

resistance.MA <- communityStab1%>%
  group_by(caseID) %>%
  filter(!'Start' %in% Stab.metric) %>%
  arrange(caseID, RD) %>%
  group_by(caseID) %>%
  slice(1) %>%
  bind_rows(., resist1) %>%
  rename(Resist.MA = totRR) %>%
  select(caseID, resp.cat, Resist.MA)

#recovery
recov.MA <- communityStab1 %>%
  filter(RD == 1) %>%
  rename(Recov.MA = totRR) %>%
  distinct(caseID, resp.cat,Recov.MA)

com.stab.MA.all <- com.stab.MA %>%
  distinct(caseID, resp.cat, AUC.delatbm.tot.MA,OEV.MA,CV.MA) %>%
  left_join(., recov.MA) %>%
  left_join(., resistance.MA)


summary(com.stab.MA.all)

#remove duplicates 
CommunityStab.MA <- distinct(com.stab.MA.all, caseID, resp.cat, OEV.MA,AUC.delatbm.tot.MA, Recov.MA, Resist.MA, CV.MA)
#rm(com.stab.MA)



####  Explore community stability ####
## see how many studies are lacking a reported community stability and why 
## e.g. response category, study design not identical

str(data.plot)
test.data <- left_join(data.plot, CommunityStab.MA, by = c( 'caseID', 'resp.cat'))

test.data %>%
  group_by(caseID, AUC.delatbm.tot.MA) %>%
  mutate(meanAUC = mean(AUC.RR, na.rm = T))%>%
ggscatter(., x = 'AUC.delatbm.tot.MA', y= 'meanAUC', cor.coef = T)

test.data <- test.data %>%
  ungroup()%>%
  drop_na(AUC.delatbm.tot.MA) %>%
  distinct(resp.cat, AUC.delatbm.tot.MA, caseID)

#### Merge community stabilit(ies) ####
# use community stability from Meta-Analysis; alternative: sum of individual species if community stability is not given 
MergedComStab <- left_join(com.stab.sum,test.data,  by = c('caseID','resp.cat')) 

write.csv(MergedComStab, file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/MergedCommunityStability.csv'))
