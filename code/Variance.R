#### Variance ####
## lnCVR over time requires a mean over time --> mean of mean is not great!!
## instead use AUC 

allData%>%
  mutate(Dist.SD = as.numeric(Dist.SD), 
         Con.SD = as.numeric(Con.SD)) %>%
  ggscatter(., x = 'Dist.M', y= 'Dist.SD', add= 'reg.line', cor.method = 'spearman', cor.coef = TRUE)

variance <- allData%>%
  distinct(caseID, studyID,  species, resp.cat, RD,Con.M, Dist.M,Con.N, Dist.N, Dist.SD, Con.SD)%>%
  mutate(dummyRR = Con.M + Dist.M) %>% 
  filter(dummyRR != 0) %>%##take out those rows where biomass is 0 in both treatment (Biomass) + control (con.bio)
  mutate(CV.Dist= (Dist.M/as.numeric(Dist.SD)),
         CV.con= (Con.M/as.numeric(Con.SD)),
         Dist.N = Dist.N, 
         Con.N = Con.N)%>%
  mutate(lnCVR= (log(CV.Dist/CV.con) + (1/(2*(Dist.N-1)))-(1/(2*(Dist.N-1)))),
         lnCV = NA) %>% #lnCV =  ln sd  ln mean 1 +(1/(2n-1))
  drop_na(lnCVR)%>%
  filter(!is.infinite(lnCVR))%>%
  group_by(caseID, species) %>%
  filter(n()>2) %>%
  summarise(AUC.lnCVR= auc(x = RD, y=lnCVR, absolutearea= F))


#### AUC Loop:Community Stability ####
#create empty df
com.stab <- tibble()

#remove duplicates from dataframe
com.response <-response %>%
  distinct(caseID, studyID, system, lat, long, organism, duration, differentiation,dist.cat,open, resp, resp.cat, RD, deltabm.tot, LRR) %>%
  mutate(Stab.metric = paste(ifelse(RD == 1, 'Recovery', ifelse(RD==0, 'Start', ''))))



