#### R script for additional analysis exploring species contributions to stability (Kunze et al. 2025 in Ecological Monographs)#

#packages
library(tidyverse)
library(readxl)
library(here)
library(metafor)
library(cowplot)
library(GGally)
library(ggpubr)
library(ggpmisc)


#### import data ####

# species stability
SpeciesStab <- read_csv('Data/SpeciesStabilities.csv')

names(SpeciesStab)
unique(SpeciesStab$organism)

#### MA - species contributions ####

#### AUC plot ####
raw <-ggplot(SpeciesStab, aes(x=AUC.pi,y=AUC.RR, col=system))+
  geom_hline(yintercept=0, colour="grey")+
  geom_vline(xintercept=0, colour="grey")+
  geom_point(alpha=0.8, size = 2)+
  scale_shape_manual(values=c(16, 17,1,2,5))+
  theme_bw()+
  labs(x = 'Relative Contribution to Stability', y = 'Absolute Contribution to Stability', color = 'System')+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(legend.title=element_text(size=13), 
        legend.text=element_text(size=12))+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(colour="black",linewidth =1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",linewidth=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.1,0.9,0.1,0.1),"cm"))
raw
ggsave(plot = raw, file = here('output/AUC.RRdeltaPi.tiff'), width = 6, height = 6)



#### Unweighted MA exploring species contributions ####
# unweighted MA requires column 1
SpeciesStab$unweighted<-1

M0<-rma.mv(AUC.RR,unweighted,
           mods = ~organism+system+dist.cat+open+resp.cat+abs(lat),
           random = ~ 1 | caseID,
           method="REML",data=SpeciesStab)

summary(M0) 

model.output <- data.frame(estimate = M0$b, # This is our vector of two effect sizes for our 2 sub-categories
                           ci.lb = M0$ci.lb, # Confidence interval lower bounds
                           ci.ub = M0$ci.ub,
                           pval = M0$pval) %>%# Confidence interval upper bounds)
rownames_to_column() %>%
rename(mods = rowname) %>%
mutate(stability = paste('Absolute Contribution'))%>%
mutate(moderators = str_replace(mods, 'organism', 'organism:'),
         moderators = str_replace(moderators, 'system', 'system:'),
         moderators = str_replace(moderators, 'dist.cat', 'dist:'),
         moderators = str_replace(moderators, 'open', 'open:'),
         moderators = str_replace(moderators, 'resp.cat', 'resp:')) 
  
model.output$moderators[model.output$moderators == 'intrcpt']<-'intercept'


##plot results
A <- ggplot(model.output, aes(x = estimate, y = moderators, color = moderators))+
  geom_vline(xintercept = 0)+
  labs(y = 'Factor', x = 'Estimate', title = 'Absolute Contribution to Stability')+
  scale_color_manual(values = c('black','black','black','darkred','black','black','black','black','black','black','black','black'))+
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), height = .2)+
  geom_point(size = 3)+
  theme_bw()+
  theme(legend.position = 'none')+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),
        axis.text.y=element_text(size=12,colour="black",angle=0,hjust=0.4),
        axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),
        axis.text.x=element_text(size=12,colour="black"),
        panel.border=element_rect(colour="black",linewidth=1.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
A

M1<-rma.mv(AUC.pi,unweighted,
           mods = ~organism+system+dist.cat+open+resp.cat+abs(lat),
           random = ~ 1 | caseID,
           method="REML",data=sppMA)

summary(M1)

model.output.pi <- data.frame(estimate = M1$b, # This is our vector of two effect sizes for our 2 sub-categories
                           ci.lb = M1$ci.lb, # Confidence interval lower bounds
                           ci.ub = M1$ci.ub,
                           pval = M1$pval) %>%# Confidence interval upper bounds)
  rownames_to_column() %>%
  rename(mods = rowname) %>%
  mutate(stability = paste('Absolute Contribution'))%>%
  mutate(moderators = str_replace(mods, 'organism', 'organism:'),
         moderators = str_replace(moderators, 'system', 'system:'),
         moderators = str_replace(moderators, 'dist.cat', 'dist:'),
         moderators = str_replace(moderators, 'open', 'open:'),
         moderators = str_replace(moderators, 'resp.cat', 'resp:')) 

model.output.pi$moderators[model.output.pi$moderators == 'intrcpt']<-'intercept'

B<- ggplot(model.output.pi, aes(x = estimate, y = moderators))+
  geom_vline(xintercept = 0)+
  labs(y = 'Factor', x = 'Estimate', title = 'Relative Contribution to Stability')+
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), height = .2)+
  geom_point(size = 3)+
  theme_bw()+
  theme(legend.position = 'none')+
  theme(axis.title.y=element_text(size=14, face="plain", colour="black",vjust=0.3),
        axis.text.y=element_text(size=12,colour="black",angle=0,hjust=0.4),
        axis.title.x=element_text(size=14,face="plain",colour="black",vjust=0),
        axis.text.x=element_text(size=12,colour="black"),
        panel.border=element_rect(colour="black",linewidth=1.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

plot_grid(A, B, labels = c('(a)', '(b)'), ncol = 2)
ggsave(plot = last_plot(), file = here('output/Significancetest.tiff'), width = 12, height = 5)


### forest ###
forest_plot <- forest(M0$b, # This is our vector of two effect sizes for our 2 sub-categories
                      ci.lb = M0$ci.lb, # Confidence interval lower bounds
                      ci.ub = M0$ci.ub, # Confidence interval upper bounds
                      annotate = TRUE, # This tells function to list effect size and CIs for each group on our graph
                      xlab = "estimate", # label for x-axis
                      slab = c('intrcpt','organism:macroalgae','organism:periphyton','organism:plant','organism:vertebrate','organism:zooplankton', 'system:marine', 'system:terrest', 'dist:phys', 'open:yes', 'duration'), #label for y-axis
                      cex = 1, # Font size for entire graph (excluding headers)
                      digits = 2 # Round effect size and CI to 2 digits
                      
)
op <- par(cex=1, font=2) # Set up font for rest of graph (just the headers of the graph remain), to make bold headings, set font=2
text(1.8, 13, "Estimate [95% CI]")
text(-2.8,13, "Mods")


