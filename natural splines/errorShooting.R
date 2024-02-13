x = c(0:25)
sinus = data.frame(x = c(0:25),y = 1*sin(0.375*x) )
ggplot(sinus, aes(x = RD, y = y))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_line()+
  theme_bw()
ggsave(plot=last_plot(), file = here('~/Desktop/phD/Meta_Multistab/Output/sinusCurveRD.png'))

sinus$RD = (sinus$x)/25

sinusLin<-myauc(sinus$RD, sinus$y,  from = min(sinus$RD, na.rm = TRUE), to = max(sinus$RD, na.rm = TRUE),
              type = c('linear'), absolutearea = FALSE)

sinusN<-myauc(sinus$RD, sinus$y,  from = min(sinus$RD, na.rm = TRUE), to = max(sinus$RD, na.rm = TRUE),
                type = c('spline'), absolutearea = FALSE)

sinusLin[[1]]
sinusN[[1]]


# Define the function to be integrated
fun <- function(x) 1 * sin(0.375 * x)

# Perform numerical integration
result <- integrate(fun, lower = 0, upper = 25)

# Extract the value of the integral (area under the curve)
area_under_curve <- result$value

# Print the area
print(area_under_curve)

#############################
hist(response$RR)
hist(stab.auc$AUC.RR.spline)

b <- filter(stab.auc, AUC.RR.spline < -2 &  AUC.RR < -0.9)

##CK045_6_Asplancha
test <- filter(response, caseID == 'CK045_4' & species == 'Cladocera')

AUC.RR<-myauc(test$DAY, test$RR,  from = min(test$DAY, na.rm = TRUE), to = max(test$DAY, na.rm = TRUE),
              type = c('linear'), absolutearea = FALSE)
AUC.pi<-myauc(test$DAY, test$delta.pi, from = min(test$DAY, na.rm = TRUE), to = max(test$DAY, na.rm = TRUE),
              type = c('linear'), absolutearea = FALSE)
AUC.RR.spline<-myauc(test$RD, test$RR,  from = min(test$RD, na.rm = TRUE), to = max(test$RD, na.rm = TRUE),
                     type = c("linear"),absolutearea = FALSE)
AUC.pi.spline<-myauc(test$RD, test$delta.pi, from = min(test$RD, na.rm = TRUE), to = max(test$RD, na.rm = TRUE),
                     type = c("linear"),absolutearea = FALSE)

a <- tibble(
       AUC.RRFit = list(AUC.RR[[2]]),
       AUC.piFit = list(AUC.pi[[2]]),
       AUC.RR = AUC.RR[[1]],
       AUC.pi = AUC.pi[[1]],
       AUC.RR.spline.Fit = list(AUC.RR.spline[[2]]),
       AUC.pi.splineFit = list(AUC.pi.spline[[2]]),
       AUC.RR.spline = AUC.RR.spline[[1]],
       AUC.pi.spline = AUC.pi.spline[[1]]
       )
a


a %>% 
  unnest(AUC.RR.spline.Fit) %>%
  rename(LinFit = yPred,
         x = xRange) %>%
  #unnest(AUC.RR.spline.Fit) %>%
  ggplot(.,)+
  #geom_hline(yintercept = 0)+
  geom_point(data = tibble(x = test$RD, y = test$RR), aes(x = x, y = y)) +
  geom_line(linetype = 'solid', aes(x = x, y = LinFit) )+
 # geom_line(linetype = 'dashed', aes(x = xRange, y = yPred))+
  theme_bw()
ggsave(plot = last_plot(), file = here('~/Desktop/phD/Meta_Multistab/Output/CK045_4_Cladocera_LinRD.png'))
AUC.RR[[1]] #-0.5464704
AUC.RR.spline[[1]] #-0.412704

############################################################################
#aucLin <- myauc(x = x, y = y, type = 'linear', absolutearea = F)

aucLin[[2]] %>% 
  ggplot(., aes(x = xRange, y = yPred)) +
  geom_point(data = tibble(x = x, y = y), aes(x = x, y = y)) +
  geom_line() +
  geom_hline(yintercept = 0)
