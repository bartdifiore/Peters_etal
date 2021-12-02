library(sf)

#---------------------------------------------------------------
## Extract excretion data at the benthic quad survey locations
#---------------------------------------------------------------


source("code/layout.R") # This pulls in an object, "polysB", which is a sp object with the locations for each quad

source("code/excretion_maps.R")

quads <- st_as_sf(polysB)
quads$transect <- coords$transect
quads$quad <- coords$quad


df <- read.csv("data/derived_excretionestimates.csv") %>% 
  st_as_sf(coords = c("dist_x", "dist_y")) %>% # convert df into a spatial data set
  group_by(plot, period)

# visualize the overlay for one site

ggplot()+
  geom_raster(data = d %>% filter(plot == "CC1", period == "Closed season"), 
              aes(x = dist_x, y = dist_y, fill = z))+
  geom_point(data = lob %>% filter(plot == "CC1", period == "Closed season"), aes(dist_x, dist_y, size = excr_total), color = "red", pch = 21)+
  geom_sf(data = quads, aes())+
  viridis::scale_fill_viridis(trans = "log1p")+
  theme_classic()


# Extract the excretion data
  # This code summarized the excretion within each quad by averaging across the extrapolated values contained in each 1 m2 quadrat. 

temp <- st_contains(quads, df[df$plot == "CC1" & df$period == "Closed season"]) # This is the function that does the work. It produces a list where the number of slots is the number of different quads (40). Just need to wrangle the list into a dataframe, where each row is a different quad, at a different site, and time point along with the associated average excretion.

#test the interspection function
temp2 <- df %>% group_by(plot, period, status) %>%
  st_join(quads, join = st_intersects) %>%
  filter(plot == "CC1", period == "Closed season") %>%
  drop_na(transect, quad)

ggplot(temp2)+
  geom_sf(aes(geometry = geometry, fill =z)) # ok so it is extracting the data from under the quads!!!


extracted_data <- df %>% group_by(plot, period, status) %>%
  st_join(quads, join = st_intersects) %>%
  drop_na(transect, quad) %>%
  ungroup() %>% 
  group_by(plot, period, status, transect, quad) %>% 
  summarize(mean_nh4 = mean(z))









ext_nh4 <- function(data){
  temp_list <- st_within(data)
  out <- vector()
  for(i in 1:dim(temp_list)[1]){
    out[i] <- mean(temp_list[[i]])
  }
  out
}

ext_nh4(df[df$plot == "CC1" & df$period == "Closed season",])

quads %>% 
  group_by(transect, quad) %>%
  summarize(mean_NH4 = mean(st_contains(df %>% filter(plot == "CC1", period == "Closed season") )))








