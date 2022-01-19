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
pl <- read.csv("data/lobster_spatial.csv") %>% janitor::clean_names() %>%
  mutate(date.d = lubridate::mdy(date)) %>%
  separate(date, into = c("month", "day", "year"), sep = "[/]")


lob <- pl %>% group_by( plot, period, status, year, month, day, date.d, sample, dist_x, dist_y) %>%
  drop_na(dist_y, dist_x) %>%
  dplyr::summarize(biomass = sum(sfdm_g), 
                   excr_total = sum(excr_ind)/60) %>% # make it per minute!!! 
  mutate(temp = ifelse(period == "Closed season", "C", "O"), 
         id = paste(plot, temp, year, month, day, sep = "-"))

ggplot(lob, aes(x = dist_x, y = dist_y))+
  geom_point(aes(size = excr_total))+
  facet_wrap(~id)


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



test <- bd_fun(lob, group = "CC1-C-19-9-14")

out <- list()
id <- unique(lob$id)

for(i in 1:length(id)){
  out[[i]] <- bd_fun(lob, group = id[i])
}

formerge <- lob %>% distinct(plot, status)

d <- bind_rows(out) %>% 
  separate(id, into = c("plot", "period", "year", "month", "day"), sep = "[-]") %>%
  mutate(period = ifelse(period == "O", "Open season", "Closed season")) %>%
  left_join(formerge)

write.csv(x = d, "data/derived_excretionestimates.csv", quote = F, row.names = F)

#------------------------------------------------------------
## Plot up the maps
#------------------------------------------------------------

# Build a multipanel plot, each panel is a heat map of the mean or SD of lobster derived nitrogen

forplot <- d %>% 
  group_by(dist_x, dist_y, plot, status) %>%
  summarize(Magnitude = mean(z), 
            sd = sd(z), 
            Consistency = 1/(sd(z)/mean(z)))

p1 <- forplot %>% filter(status == "Reserve") %>%
  pivot_longer(cols = c(Magnitude, Consistency )) %>%
  mutate(name = forcats::fct_rev(name)) %>%
  ggplot(aes(x = dist_x, y = dist_y))+ # plot it up! 
  geom_raster(aes(fill = value))+
  viridis::scale_fill_viridis(trans = "log1p")+
  facet_grid(plot ~ name)+
  labs(x = "", y = "", title = "Reserve sites")+
  theme_classic()

ggsave("figures/reserve.png", p1, width = 8.5, height = 8.5)

p2 <- forplot %>% filter(status == "Off reserve") %>%
  pivot_longer(cols = c(Magnitude, Consistency)) %>%
  mutate(name = forcats::fct_rev(name)) %>%
  ggplot(aes(x = dist_x, y = dist_y))+ # plot it up! 
  geom_raster(aes(fill = value))+
  viridis::scale_fill_viridis(trans = "log1p")+
  facet_grid(plot ~ name)+
  labs(x = "", y = "", title = "Off reserve sites")+
  theme_classic()
ggsave("figures/offreserve.png", p2, width = 8.5, height = 8.5)

# Build a multipanel plot of magnitude for each plot before and after season opened

p3 <- d %>%
  group_by(dist_x, dist_y, plot, period, status) %>%
  summarize(Magnitude = mean(z)) %>%
  ggplot(aes(dist_x, dist_y))+
  geom_raster(aes(fill = Magnitude))+
  viridis::scale_fill_viridis(trans = "log1p")+
  facet_grid(plot ~ period)+
  labs(x = "", y = "")+
  theme_classic()

ggsave("figures/magnitudeXperiod.png", p3, width = 8.5, height = 2*8.5)





ggplot(d, aes(x = dist_x, y = dist_y))+ # plot it up! 
  geom_raster(aes(fill = z))+
  geom_point(data = lob, aes(dist_x, dist_y, size = excr_total), color = "red", pch = 21)+
  viridis::scale_fill_viridis(trans = "log1p")+
  facet_wrap(~id)+
  theme_classic()

ggsave("figures/lobster_maps.png", width = 15, height = 12)


id <- unique(lob$id)
for(i in id){
  temp_plot <- ggplot(data = subset(d, plot == i), aes(x = dist_x, y = dist_y))+ # plot it up! 
    geom_raster(aes(fill = z))+
    geom_point(data = subset(lob, plot == i), aes(dist_x, dist_y, size = biomass), color = "red", pch = 21)+
    viridis::scale_fill_viridis(trans = "log1p")+
    labs(title = paste0("plot_", i))+
    theme_classic()
  
  ggsave(temp_plot, path = "figures/", file=paste0("plot_", i,".png"), width = 14, height = 10, units = "cm")
  
}











