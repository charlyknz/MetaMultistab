### correlation ###
AllStab_wide <-  AllStab %>%
  filter(!str_detect(RD.metric,'abs' )) %>%
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
#ggsave(plot = last_plot(), file = here('output/Correlation_RealisedResponseTraits_Instab.png'), width = 8, height = 7)


##### RD and MA Community Stability#####
str(ComStab)
AllStab1 <- merge(realised.pert,MergedComStab, by = c('caseID', 'resp.cat')) %>%
  filter(!str_detect(RD.metric,'abs' ))%>%
  select(caseID, resp.cat, RD.metric, RD.value, AUC.delatbm.tot.MA ) %>%
  spread(key = RD.metric, value = RD.value)
str(AllStab1)

mean2<- ggscatter(AllStab1, x = 'mean_spp_deltabm', y = 'AUC.delatbm.tot.MA', 
                  xlab = 'Mean species response trait', ylab = '(reported) Relative OEV',
                  add = 'reg.line', cor.coef = T, cor.method='spearman')
mean2

RD_div2<- ggscatter(AllStab1, x = 'RD_div', y = 'AUC.delatbm.tot.MA', 
                    xlab = 'Response divergence', ylab = '(reported) Relative OEV',
                    add = 'reg.line', cor.coef = T, cor.method='spearman')
RD_div2

RD_diss2<- ggscatter(AllStab1, x = 'RD_diss', y = 'AUC.delatbm.tot.MA', 
                     xlab = 'Response dissimilarity', ylab = '(reported) Relative OEV',
                     add = 'reg.line', cor.coef = T, cor.method='spearman')
RD_diss2

plot_grid(RD_div2, RD_diss2,mean2, labels = c('(a)', '(b)', '(c)'))
#ggsave(plot = last_plot(), file = here('output/Correlation_RealisedResponseTraits_reportedInstab.png'), width = 8, height = 7)

