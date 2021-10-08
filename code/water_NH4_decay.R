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


#--------------------------------------------------------------------------------
## Excretion maps
## Bart DiFiore
## Code originally written: October, 6 2021
#--------------------------------------------------------------------------------


# Build maps based on Joey's model

# So Joey's model essentially says that the y.intercept of the model (parameter d) if determined by the abundance of lobsters. The abundance of lobsters is a linear predictor of the excretion rates. Therefore, we use the slope (parameter e) to predict the decay relationship (e.g. ~ distance) based on the idea that the total amount of ammonia produced in one hour is the concentration at zero. There was no difference in slope with the abundance of lobster, which simplifies matters. We will just use the total ammonia excreted in one hour as the concentration at zero and model the decay relationship with distance. There as also no difference with lobster abundance in the horizonal asymptote of the function. So similarly, I've used the mean across abudnance groups for the slope (e) and horizonal asymptote (c) parameters. To make sure that we are reaching zero, I will subtract the mean horizonal asymptote value from all predictions, so that we are correcting for background nitrogen production. 

fun <- function(distance, value.0, slope = 0.6302233, h.asymptote = 0.93361){
  raw  = h.asymptote + (value.0-h.asymptote)*(exp(-distance/slope))
  raw - h.asymptote # correct for background nitrogen
}

# Example
df <- expand.grid(distance = seq(1, 4, length.out = 100), value.0 = c(100, 200, 400))
df$y <- fun(distance = df$distance, value.0 = df$value.0)

ggplot(df, aes(x = distance, y = y))+
  geom_line(aes(color = as.factor(value.0)))


# Build the maps
pl <- read.csv("data/lobster_spatial.csv") %>% janitor::clean_names()


lob <- pl %>% group_by(period, status, date, sample, plot, dist_x, dist_y) %>%
  drop_na(dist_y, dist_x) %>%
  dplyr::summarize(biomass = sum(sfdm_g), 
                   excr_total = sum(excr_ind)/60) %>% # make it per minute!!! 
  mutate(temp = ifelse(period == "Closed season", "C", "O"), 
         id = paste(plot, temp, sep = "-"))

ggplot(lob, aes(x = dist_x, y = dist_y))+
  geom_point(aes(size = excr_total))+
  facet_wrap(~plot+period)


# Build a function to estimate nitrogen at every point in the grid. 

grid <- expand.grid(dist_x = seq(0, 20, length.out = 100), dist_y = seq(0,20, length.out = 100)) # construct a grid of points across which to estiamte nitrogen
plot(grid) # just confirm

lob <- as.data.frame(lob)


bd_fun <- function(data,group){
  dat <- data[data$id == group,]
  points <- dat[,c("dist_x", "dist_y")]
  excr_total <- dat[, c("excr_total")]
  out <- t(fields::rdist(points, grid)) # estimate the euclidean distance from each focal point to each other point in the grid, transpose such that the matrix is organized where each row is a different grid point, and each column is a different focal point
  
  temp <- matrix(nrow = dim(out)[1], ncol = dim(out)[2]) # build a matrix to store the output of the for loop
  
  
  for(i in 1:dim(out)[2]){
    temp[,i] <- fun(distance = out[,i], value.0 = excr_total[i]) # predict the contribution to the z variable from each focal point to each other point in the grid
  }
  
  z <- temp %*% matrix(data = 1, nrow = dim(out)[2], ncol = 1) # sum across rows in the matrix, such that it is the sum of the contributions to Z from each focal point where lobsters are located
  
  grid$z <- z # return Z
  grid$id <- group
  grid
  
}



test <- bd_fun(lob, group = "CC1-C")

out <- list()
id <- unique(lob$id)

for(i in 1:length(id)){
  out[[i]] <- bd_fun(lob, group = id[i])
}


d <- bind_rows(out) %>% 
  left_join(distinct(dplyr::select(lob, c(period, status, id)))) %>% 
  separate(id, into = c("plot", "junk"), sep = "[-]")
  
  
  
  res <- df %>% mutate(category=cut(a, breaks=c(-Inf, 0.5, 0.6, Inf), labels=c("low","middle","high")))

d %>% group_by(plot, period) %>% 
  summarize(mean = mean(z, na.rm = T))


ggplot(d, aes(x = dist_x, y = dist_y))+ # plot it up! 
  geom_raster(aes(fill = z))+
  geom_point(data = lob, aes(dist_x, dist_y, size = excr_total), color = "red", pch = 21)+
  viridis::scale_fill_viridis(trans = "log1p")+
  facet_wrap(~plot+period)+
  theme_classic()

ggsave("figures/lobster_maps.png", width = 15, height = 12)


plot <- unique(lob$plot)
for(i in plot){
  temp_plot <- ggplot(data = subset(d, plot == i), aes(x = dist_x, y = dist_y))+ # plot it up! 
    geom_raster(aes(fill = z))+
    geom_point(data = subset(lob, plot == i), aes(dist_x, dist_y, size = biomass), color = "red", pch = 21)+
    viridis::scale_fill_viridis(trans = "log1p")+
    labs(title = paste0("plot_", i))+
    facet_wrap(~period)+
    theme_classic()
  
  ggsave(temp_plot, path = "figures/", file=paste0("plot_", i,".png"), width = 14, height = 10, units = "cm")
  
}

temp <- d %>% dplyr::select(-junk) %>% group_by(plot, period, dist_x, dist_y) %>%
  pivot_wider(names_from = period, values_from = z) %>%
  mutate(diff.z = (`Open season` - `Closed season`) / (`Closed season` + `Open season`))


ggplot(temp, aes(x = dist_x, y = dist_y ))+
  geom_raster(aes(fill = diff.z))+
  #scale_fill_gradient2( low = "blue", mid = "white", high = "red")+
  scale_fill_distiller(type = "div", palette = "RdBu", direction = 1)+
  facet_wrap(~plot+status)+
  theme_classic()

ggsave(filename = "figures/difference.png", width = 10, height = 8)



# Statistical analyses


temp <- d %>% dplyr::select(-junk) %>% group_by(plot, period, dist_x, dist_y) %>%
  pivot_wider(names_from = period, values_from = z) %>%
  mutate(diff.z = (`Open season` - `Closed season`) / (`Closed season` + `Open season`)) %>%
  rename(open = `Open season`, closed = `Closed season`)

ggplot(temp, aes(x = closed, y = open))+
  geom_point(aes(color = plot, shape = status))+
  facet_wrap(~ plot)

ggplot(d, aes(x = period, y = z))+
  geom_boxplot(aes(color = status))

library(lme4)
library(lmerTest)
mod.lmer <- lmer(z ~ period*status + (1|plot), d)
summary(mod.lmer)
car::qqPlot(residuals(mod.lmer))
hist(residuals(mod.lmer))


ls_means(mod.lmer, which = c(""))
emmeans::emmeans(mod.lmer,  ~ period * status,)


pred <- ggeffects::ggpredict(mod.lmer, terms = ~period*status)

plot(pred)


mod.glmer <- glmer(z ~ period*status + (1|plot), d, family = Gamma(link = "log"))
summary(mod.glmer)
car::qqPlot(residuals(mod.glmer))
hist(residuals(mod.glmer))
pred <- ggeffects::ggpredict(mod.glmer, terms = ~period*status)
plot(pred)






