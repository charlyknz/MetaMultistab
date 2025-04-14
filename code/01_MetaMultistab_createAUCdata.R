#### R script for meta analysis ####

library(tidyverse)
library(readxl)
library(MESS)
library(here)
library(cowplot)
library(GGally)
library(ggpubr)

#### import data ####
study <- read_excel("~/Desktop/phD/Meta_Multistab/MetaMultistab/Data/Multistab_species_data.xlsx") %>%
  select(-c(35:51)) 
names(study)

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

#remove duplicates (if there are any)
CommunityStab.MA <- distinct(com.stab.MA.all, caseID, resp.cat, OEV.MA,AUC.delatbm.tot.MA, Recov.MA, Resist.MA, CV.MA)

#write_csv(CommunityStab.MA, file = here('~/Desktop/phD/Meta_Multistab/MetaMultistab/output/CommunityStabilities.csv'))
  

