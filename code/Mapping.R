#--------------------------------------------------------------------------------
## Excretion maps
## Bart DiFiore
## Code originally written: October, 6 2021
#--------------------------------------------------------------------------------

# The goal here is to estimate a statistical function describing the decline in nitrogen concentration as a function of lobster abundance and distance from den sites. Using this regression relationship, this code will then estimate the nitrogen concentration as a function of geographic space and lobster abundance. 

#---------------------------------------------------------------------------------
## Setup
#---------------------------------------------------------------------------------

library(tidyverse)

#---------------------------------------------------------------------------------
## Experimental regression relationship
#---------------------------------------------------------------------------------

df <- read.csv("data/water_nutrients.csv") %>% janitor::clean_names()

ggplot(df, aes(x = distance, y = nh4_u_m))+
  geom_point(aes(size = count, color = site))


lm1 <- lm(log(nh4_u_m) ~ log(distance+0.001) * log(sfdm_g), df)
summary(lm1)

newdat <- expand.grid(distance = seq(0,4,length.out = 100), sfdm_g = mean(df$sfdm_g))
newdat$nh4_u_m <- exp(predict(lm1, newdata = newdat))

ggplot(df, aes(x = distance+0.001, y = nh4_u_m))+
  geom_point(aes(size = sfdm_g, color = site))+
  geom_line(data = newdat, aes(x = distance, y = nh4_u_m, group = sfdm_g))+
  scale_y_log10()+
  scale_x_log10()

ggplot(df, aes(x = distance+0.001, y = nh4_u_m))+
  geom_point(aes(size = sfdm_g, color = site))+
  geom_line(data = newdat, aes(x = distance, y = nh4_u_m, group = sfdm_g))

ggplot(df, aes(x = as.factor(distance), y = nh4_u_m))+
  geom_boxplot()

# Ok so these models suggest that nitrogen concentration does not approach zero. This isn't biologically realistic, as nitrogen from the lobster has to approach zero with distance. The measures are picking up background nitrogen in the water column. Therefore, I am going to offset the measures by the average nitrogen concentration at sites with one lobster, at the farthest distance. Therefore, we assume that the contributions to nitrogen from lobster are zero at 4 m, if only 1 lobter is present. To keep values from going negative, I assumed if the offset is negative, ammonium is zero. This may or may not be realistic. An alternative would be to use the lower asymptote of the function as the offset. 

df %>%
  filter(count == 1) %>%
  ggplot(aes(x = distance, y = nh4_u_m))+
  geom_point(aes(color = site))

avg <- df %>%
  filter(count == 1, distance == 4) %>%
  summarize(avg = mean(nh4_u_m))
avg <- as.vector(avg$avg)

df <- df %>%
  mutate(nh4_off = ifelse(nh4_u_m - avg > 0, nh4_u_m -avg, 0))


lm2 <- lm(log(nh4_off+1) ~ log(distance+1) * log(sfdm_g), df)
summary(lm2)
plot(residuals(lm2) ~ log(distance), df)

newdat <- expand.grid(distance = seq(0,4,length.out = 100), sfdm_g = mean(df$sfdm_g))
newdat$nh4_off <- exp(predict(lm2, newdata = newdat))

ggplot(df, aes(x = log(distance+1), y = log(nh4_off+1)))+
  geom_point(aes(size = sfdm_g, color = site))+
  geom_line(data = newdat, aes(x = log(distance+1), y = log(nh4_off+1)))
  # scale_y_log10()+
  # scale_x_log10()

ggplot(df, aes(x = distance, y = nh4_off))+
  geom_point(aes(size = sfdm_g, color = site))+
  geom_line(data = newdat, aes(x = distance, y = nh4_off, group = sfdm_g))


# So this model feels adequate for now (R2 ~ 0.75). Certainly not perfect and we could refine this model, so that the initial decline isn't so steep using some other nonlinear functions. 

#----------------------------------------------------------------------------------
## Extrapolate to one site
#----------------------------------------------------------------------------------

pl <- read.csv("data/lobster_spatial.csv") %>% janitor::clean_names()


lob <- pl %>% group_by(period, status, date, sample, plot, dist_x, dist_y) %>%
  drop_na(dist_y, dist_x) %>%
  dplyr::summarize(biomass = sum(sfdm_g), 
                   excr_total = sum(excr_ind)) %>% 
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
  biomass <- dat[, c("biomass")]
  out <- t(fields::rdist(points, grid)) # estimate the euclidean distance from each focal point to each other point in the grid, transpose such that the matrix is organized where each row is a different grid point, and each column is a different focal point
  
  temp <- matrix(nrow = dim(out)[1], ncol = dim(out)[2]) # build a matrix to store the output of the for loop
  
  
  for(i in 1:dim(out)[2]){
    temp[,i] <- exp(predict(lm2, newdata = data.frame(distance = out[,i], sfdm_g = biomass[i]))) # predict the contribution to the z variable from each focal point to each other point in the grid
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
  left_join(distinct(select(lob, c(period, status, id)))) %>% 
  separate(id, into = c("plot", "junk"), sep = "[-]")

d %>% group_by(id) %>% 
  summarize(mean = mean(z, na.rm = T))


ggplot(d, aes(x = dist_x, y = dist_y))+ # plot it up! 
  geom_raster(aes(fill = z))+
  geom_point(data = lob, aes(dist_x, dist_y, size = biomass), color = "red", pch = 21)+
  viridis::scale_fill_viridis()+
  facet_wrap(~id+status)+
  theme_classic()


plot <- unique(lob$plot)
for(i in plot){
  temp_plot <- ggplot(data = subset(d, plot == i), aes(x = dist_x, y = dist_y))+ # plot it up! 
    geom_raster(aes(fill = z))+
    geom_point(data = subset(lob, plot == i), aes(dist_x, dist_y, size = biomass), color = "red", pch = 21)+
    viridis::scale_fill_viridis(limits = c(0, 100))+
    labs(title = paste0("plot_", i))+
    facet_wrap(~period)+
    theme_classic()
  
  ggsave(temp_plot, file=paste0("plot_", i,".png"), width = 14, height = 10, units = "cm")
  
}

temp <- d %>% select(-junk) %>% group_by(plot, period, dist_x, dist_y) %>%
  pivot_wider(names_from = period, values_from = z) %>%
  mutate(diff.z = (`Closed season` - `Open season`) / `Closed season` )


ggplot(temp, aes(x = dist_x, y = dist_y ))+
  geom_raster(aes(fill = diff.z))+
  scale_fill_gradient( low = "gray", high = "red")+
  facet_wrap(~plot+status)+
  theme_classic()
  











