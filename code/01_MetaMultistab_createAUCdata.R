#### R script: data wrangling for meta analysis ####

library(tidyverse)
library(readxl)
library(MESS)
library(here)
library(cowplot)
library(GGally)
library(ggpubr)

### to start the analysis please download the meta-analysis data on figshare and store in data folder###
dir.create(here('Data')) 

### store output ###
dir.create(here('output')) 

#### import data ####
study <- read_excel("Data/Multistab_species_data.xlsx") %>%
  select(-c(35)) 
names(study)

rawData <- read_excel("Data/Multistab_species_data.xlsx", 
                      sheet = "species data")

allData <- study%>%
  select(-func, -resp.cat,-resp) %>%
  merge(., rawData,by = c('caseID'))%>%
  filter(spec.inf %in% c('species', 'taxa')) %>%
  drop_na(caseID)

which(is.na(allData$Con.M))


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


#### AUC Loop species stability ####
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

### check output ###
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

write_csv(data.plot, file = here('Data/SpeciesStabilities.csv'))



#### AUC Loop Community Stability data ####
names(response)
unique(response$caseID)

communityStab1<- response %>%
  distinct(caseID,duration,DAY,RD,resp.cat,deltabm.tot)%>%
  drop_na(deltabm.tot) 

communityStab1$resp.cat[communityStab1$resp.cat == 'cover'] <-'biomass'

## create USI to run loop
USIc <- unique(communityStab1$caseID)
names(communityStab1)

#empty df
com.stab <- data.frame()

for(i in 1:length(USIc)){
  temp<-communityStab1[communityStab1$caseID==USIc[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>3){#does the next step only if at least 3 data points are present
    OEV<-auc(temp$RD, temp$deltabm.tot, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                    type = c("linear"),absolutearea = TRUE)
    AUC.delatbm.tot  <-auc(temp$RD, temp$deltabm.tot, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                               type = c("linear"),absolutearea = FALSE)
    CV<- mean(temp$deltabm.tot, na.rm = T)/sd(temp$deltabm.tot, na.rm = T) # coefficient of variation
    com.stab<-rbind(com.stab,
                      data.frame(temp[1,c(1)],
                                 OEV,CV,AUC.delatbm.tot))
    rm(temp)
  }
}
str(com.stab)
unique(com.stab$caseID)

#### Stability metrics community ####
StabMetrics <- communityStab1 %>%
  select(caseID, DAY, RD, deltabm.tot) 

### For resistance we have to slice the first entry after the start community ##
# as sometimes we dont have information on the start community, we will split the two df 
resist1 <- StabMetrics%>%
  arrange(caseID, RD) %>%
  group_by(caseID) %>%
  slice(1) %>%
  filter(RD != 0)
  
resist2 <- StabMetrics%>%
  filter(RD == 0)%>%
  distinct(caseID) %>%
  left_join(.,StabMetrics)%>%
  arrange(caseID, RD) %>%
  group_by(caseID) %>%
  slice(2)

resistance <- resist1%>%
  bind_rows(., resist2) %>%
  rename(Resist = deltabm.tot) %>%
  select(caseID, Resist)

#recovery
recov.MA <- communityStab1 %>%
  filter(RD == 1) %>%
  rename(Recov = deltabm.tot) %>%
  distinct(caseID,Recov)

com.stab.MA.all <- recov.MA %>%
  left_join(., resistance) %>%
  right_join(., com.stab) 
  
summary(com.stab.MA.all)


write_csv(com.stab.MA.all, file = here('Data/CommunityStabilities.csv'))

