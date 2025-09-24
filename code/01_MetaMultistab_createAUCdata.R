#### R script: data wrangling for meta analysis ####

library(tidyverse)
library(readxl)
library(MESS)
library(here)
library(cowplot)
library(GGally)
library(ggpubr)

### to start the analysis please download the meta-analysis data on figshare and store in data folder###
#dir.create(here('Data')) 

### store output ###
#dir.create(here('output')) 

#### import data ####
study <- read_excel("Data/Multistab_species_data_repro.xlsx", 
                    col_types = c("text", "text", "text", 
                                  "text", "text", "text", "text", "numeric", 
                                  "numeric", "text", "numeric", "text", 
                                  "numeric", "numeric", "text", "text", 
                                  "text", "text", "text", "text")) %>% 
  rename(comment = "...20")
names(study)

rawData <- read_excel("Data/Multistab_species_data_repro.xlsx", 
                      sheet = "species data", col_types = c("text", 
                                                            "text", "text", "text", "text", "text", 
                                                            "numeric", "numeric", "text", "text", 
                                                            "text", "text", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric"), 
                      na = "NA")
str(rawData)

#new studies since 2017
dataSpp <- read_excel("Data/ELE_Review_MetaMultistab_since2017.xlsx", 
                                 sheet = "species_data")

meta <- read_excel("Data/ELE_Review_MetaMultistab_since2017.xlsx", 
                      sheet = "meta") %>% 
  mutate(lat = as.numeric(lat))
unique(meta$studyID)
#merge

metadata <- study%>%
  bind_rows(., meta) %>% 
  select(-func, -resp) 

allData <- rawData %>% 
  bind_rows(., dataSpp) %>% 
  select(-comment) %>% 
  merge(., metadata, by = c("caseID", "resp.cat") )%>%
 # filter(spec.inf %in% c('species', 'taxa')) %>%
  drop_na(caseID)

which(is.na(allData$Con.M)) 
is.na(allData$Con.M)<- 0


# test if all data are imported correctly
unique(rawData$caseID)
unique(allData$studyID)
unique(allData$caseID)


setdiff( study$caseID,allData$caseID) 
# look at data
unique(allData$organism)

#remove negative biomass and abundance entries!
allData$Con.M[(allData$Con.M <0)]<-0
allData$Dist.M[(allData$Dist.M <0)]<-0

response <- allData %>%
  select(caseID, studyID, spec.inf, organism, system, lat, long, duration, differentiation, 
         species, species_specification,func,resp, resp.cat,DAY, RD,Con.M, Dist.M,Con.N, Dist.N, Dist.SD, Con.SD)%>%
  mutate(dummyRR = Con.M + Dist.M) %>% 
  filter(dummyRR != 0) %>%##take out those rows where biomass is 0 in both treatment (Biomass) + control (con.bio)
  group_by(caseID, organism, resp.cat,RD)%>%
  mutate(con.tot = sum(Con.M, na.rm = T),
         dist.tot = sum(Dist.M, na.rm = T))%>%
  ungroup() %>%
  mutate(dist.pi = Dist.M/dist.tot,
         con.pi = Con.M/con.tot,
         delta.pi = dist.pi - con.pi, # calculate effect sizes
         RR = (Dist.M-Con.M)/(Dist.M+Con.M),
         diff = (Dist.M-Con.M),
         LRR = log(dist.tot/ con.tot),
         var.lrr = ((Dist.SD)^2/Dist.N*(Dist.M)^2)+((Con.SD)^2/Con.N*(Con.M)^2),
         deltabm.tot = (dist.tot - con.tot)/(dist.tot+con.tot)) %>%
  mutate(USI = paste(caseID, resp.cat,organism, species,  sep = "_"))  %>%
  filter(resp.cat != "contribution to production") %>%
  distinct(caseID, studyID, spec.inf, organism, system, lat, long, duration, differentiation, species, species_specification,func,resp, resp.cat,DAY, RD,Con.M, Dist.M,Con.N, Dist.N,
           deltabm.tot,LRR,var.lrr,RR,diff,delta.pi,con.pi,dist.tot, con.tot,dist.pi,USI)


#remove infinite values for species relative proportion delta.pi and absolute biomass response ratio rr
response$delta.pi[is.infinite(response$delta.pi)]<-NA
response$RR[is.infinite(response$RR)]<-NA
response$RR[response$RR == 'NaN']<-NA
hist(response$RR)

#### AUC Loop species stability ####
USI <- response$USI #unique identifier

#order after time steps
response <- response[order(response$RD),]
names(response)

#create empty df
stab.auc <- tibble()

for(i in 1:length(USI)){
  temp<-response[response$USI==USI[i], ]#creates a temporary data frame for each case
   if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.RR<-auc(temp$RD, temp$RR,  from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                absolutearea = FALSE)
    AUC.diff<-sum(temp$RR, na.rm = T)
   # AUC.pi<-auc(temp$RD, temp$delta.pi, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
    #            absolutearea = FALSE)
    con.pi <- mean(temp$con.pi, na.rm = T)
   # mean.delta.pi <- mean(temp$delta.pi, na.rm = T)
    mean.RR <- mean(temp$RR, na.rm = T)
    sum.con <- sum(temp$Con.M)
     stab.auc<-rbind(stab.auc,
                    tibble(temp[1,c(1:14)],
                           AUC.RR ,
                          # AUC.pi ,
                           con.pi,
                           AUC.diff,
                           sum.con,
                           mean.RR#, mean.delta.pi
                          ))
    rm(temp)
  }
}

### check output ###
unique(stab.auc$caseID)
str(stab.auc)

data.plot <- stab.auc %>%
  distinct(caseID, studyID,system, lat, organism,duration, species, resp.cat, sum.con,AUC.RR, con.pi,mean.RR)%>%
  filter(resp.cat %in% c('abundance', 'biomass') ) %>%
  ungroup() %>%
 # drop_na(AUC.pi)%>% 
  #remove studies with less than 2 species # Note: sometimes species get removed during the AUC loop
  group_by(caseID) %>% 
  filter(n() >1)
#rm(stab.auc)

names(data.plot)
hist(data.plot$AUC.RR)

write_csv(data.plot, file = here('Data/SpeciesStabilities.csv'))



#### AUC Loop Community Stability data ####
names(response)
unique(response$caseID)

communityStab1<- response %>%
  distinct(caseID,duration,DAY,RD,resp.cat,organism,deltabm.tot)%>%
  mutate(USI = paste(caseID, resp.cat, organism, sep = "_")) 

communityStab1$resp.cat[communityStab1$resp.cat == 'cover'] <-'biomass'

## create USI to run loop
USIc <- unique(communityStab1$USI)
names(communityStab1)

#empty df
com.stab <- data.frame()

for(i in 1:length(USIc)){
  temp<-communityStab1[communityStab1$USI==USIc[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    OEV<-auc(temp$RD, temp$deltabm.tot, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                    type = c("linear"),absolutearea = TRUE)
    AUC.delatbm.tot  <-auc(temp$RD, temp$deltabm.tot, from = min(temp$RD, na.rm = TRUE), to = max(temp$RD, na.rm = TRUE),
                               type = c("linear"),absolutearea = FALSE)
    CV<- mean(temp$deltabm.tot, na.rm = T)/sd(temp$deltabm.tot, na.rm = T) # coefficient of variation
    com.stab<-rbind(com.stab,
                      data.frame(temp[1,c(1,5,6,8)],
                                 OEV,CV,AUC.delatbm.tot))
    rm(temp)
  }
}
str(com.stab)
unique(com.stab$caseID)

 #### Stability metrics community ####
StabMetrics <- communityStab1 %>%
  select(caseID, resp.cat,organism, DAY, RD, deltabm.tot) 

### For resistance we have to slice the first entry after the start community ##
# as sometimes we dont have information on the start community, we will split the two df 
resist1 <- StabMetrics%>%
  arrange(caseID, resp.cat,organism, RD) %>%
  group_by(resp.cat,organism, caseID) %>%
  slice(1) %>%
  filter(RD != 0)
  
resist2 <- StabMetrics%>%
  filter(RD == 0)%>%
  distinct(resp.cat,organism, caseID) %>%
  left_join(.,StabMetrics)%>%
  arrange(caseID, resp.cat,organism, RD) %>%
  group_by(resp.cat,organism, caseID) %>%
  slice(2)

resistance <- resist1%>%
  bind_rows(., resist2) %>%
  rename(Resist = deltabm.tot) %>%
  select(caseID, resp.cat,organism, Resist)

# resilience
resil1 <- StabMetrics%>%
  arrange(caseID,resp.cat, organism, RD) %>%
  group_by(resp.cat,organism, caseID) %>%
  slice(-c(1)) %>%
  filter(RD != 0)

resil2 <- StabMetrics%>%
  filter(RD == 0)%>%
  distinct(resp.cat,organism, caseID) %>%
  left_join(.,StabMetrics)%>%
  arrange(resp.cat,caseID,organism,  RD) %>%
  group_by(resp.cat,organism, caseID) %>%
  slice(-c(1:2))

resilience <- resil1%>%
  bind_rows(., resil2) %>%
  rename(Resil = deltabm.tot) %>%
  select(caseID,resp.cat,RD,organism,  Resil) %>% 
  mutate(USI = paste(caseID, resp.cat, organism, sep ="_"))

caseID <- unique(resilience$USI)
#create an empty data frame
slope.RD<-tibble()

# the following loop cycles through all unique cases

for(i in 1:length(caseID)){
  temp<-resilience[resilience$USI==caseID[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    lm1<-lm(Resil~log(RD+1), temp)#makes a linear regreassion
    intcp.lm.RD <- coef(summary(lm1))[1, 1]#selects the intercept
    se.intcp.lm.RD<- coef(summary(lm1))[1, 2]#selects its standard error
    resil.lm.RD <- coef(summary(lm1))[2, 1]#selects the slope
    se.slp.lm.RD<- coef(summary(lm1))[2, 2]#selects its standard error
    sd.res.lm.RD<- sd(resid(lm1)) #selects the standard deviation of the residuals
    temp.stab.lm.RD<-1/sd.res.lm.RD
    p.lm.RD<-anova(lm1)$'Pr(>F)'[1]#gives the p-value
    slope.RD<-rbind(slope.RD,data.frame(temp[1,"caseID"],intcp.lm.RD, se.intcp.lm.RD, resil.lm.RD, se.slp.lm.RD, sd.res.lm.RD,temp.stab.lm.RD,p.lm.RD))
    rm(temp)
  }
}

names(slope.RD)
summary(slope.RD)


#recovery
recov.MA <- communityStab1 %>%
  filter(RD == 1) %>%
  rename(Recov = deltabm.tot) %>%
  distinct(caseID,resp.cat,organism, Recov)

com.stab.MA.all <- recov.MA %>%
  left_join(., resistance) %>%
  right_join(., com.stab)  %>% 
  left_join(., slope.RD) %>% 
  select(-c(se.intcp.lm.RD, se.slp.lm.RD,sd.res.lm.RD, p.lm.RD, intcp.lm.RD))
  
summary(com.stab.MA.all)


write_csv(com.stab.MA.all, file = here('Data/CommunityStabilities.csv'))

