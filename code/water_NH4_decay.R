library(tidyverse)


####WATER NUTRIENTS#####
water_chem <- read.csv("data/water_nutrients.csv",stringsAsFactors = F,na.strings=".")
water_chem$Density [water_chem$Count %in% c(12,16,20)] <-"High"
water_chem$Density [water_chem$Count %in% c(6,8,9)] <-"Medium"
water_chem$Density [water_chem$Count %in% c(1,3,4)] <-"Low"


lm1 <- lm(log(NH4_uM) ~ log(Distance+0.001)*Density, water_chem)
summary(lm1)

library(drc)
model1 <- drm(NH4_uM ~ Distance, Density, fct=EXD.3(names=c("h.asymtope", "y.intercept", "Slope")),
              data = water_chem)
summary(model1)
#compare parameters
compParm(model1, "h.asymtope") # Ok so they are all getting to the same place (i.e. they are asymptotying at the background nitrogen level)
compParm(model1, "y.intercept") # How much at zero is dependent on how many lobsters... exactly what we would expect
compParm(model1, "Slope") # And they all decline at the same rate. So how much you start with doesn't change how fast it declines meaning that we could use a single slope estimate.
modelFit(model1)

#assumptions
qqnorm(resid(model1))
plot(resid(model1)~Distance, water_chem)


newdat <- expand.grid(Distance = seq(0,4,length.out = 100), Density = water_chem$Density)
newdat$NH4_uM <- predict(model1, newdata = newdat)

#reorder factor levels for plot
water_chem$Density <- factor(water_chem$Density,
                             levels = c('High', 
                                        'Medium', 
                                        'Low'))
newdat$Density <- factor(newdat$Density,
                         levels = c('High', 
                                    'Medium', 
                                    'Low'))

ggplot(water_chem, aes(x = Distance, y = NH4_uM))+
  geom_point(aes(fill = Density), size=4, pch=21)+
  geom_line(data = newdat, aes(x = Distance, y = NH4_uM, colour = Density), lwd=1.2)+
  labs(x = "Distance from lobster den (m)", 
       y = bquote(mu*M ~NH[4]^`+`))+  
  # viridis::scale_fill_viridis(option="plasma")+
  #scale_fill_gradient2(midpoint = 1,
  #                    low = "red",
  #                   mid = "yellow",
  #                  high = "blue") +
  scale_fill_manual(values = c("#e34a33","#fec44f", "#9ebcda")) +
  scale_colour_manual(values = c("#e34a33","#fec44f", "#9ebcda")) +
  theme_classic() +
  theme(axis.text.y = element_text(family = "Helvetica", size=12, colour='black'),
        axis.title.y = element_text(family = "Helvetica", size=14,  colour='black'), 
        axis.text.x =element_text(family = "Helvetica", size=12,  colour='black'),
        axis.title.x = element_text(family = "Helvetica", size=14,  colour='black'),
        legend.title = element_text(family = "Helvetica", size=12,  colour='black', face="bold"),
        legend.text = element_text(family = "Helvetica", size=12,  colour='black'),
        aspect.ratio=1)

ggsave("water_nutrients.pdf", width=5, height=5, dpi=300)

#mean concentrations by dens density and distance
water_chem_summary <- water_chem %>%
  dplyr:: group_by (Distance, Density) %>% 
  dplyr:: summarise(N = length(Sample_ID), 
                    MEAN = mean(NH4_uM, na.rm=TRUE),
                    SE=sd(NH4_uM)/sqrt(N))

#mean concentrations by distance
water_chem_summaryb <- water_chem %>%
  dplyr:: group_by (Distance) %>% 
  dplyr:: summarise(N = length(Sample_ID), 
                    MEAN = mean(NH4_uM, na.rm=TRUE),
                    SE=sd(NH4_uM)/sqrt(N))
