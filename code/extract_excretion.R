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
  geom_raster(data = d %>% filter(plot == "CC1", period == "Closed season", sample == "A"), 
              aes(x = dist_x, y = dist_y, fill = z))+
  geom_point(data = lob %>% filter(plot == "CC1", period == "Closed season"), aes(dist_x, dist_y, size = excr_total), color = "red", pch = 21)+
  geom_sf(data = quads, aes())+
  viridis::scale_fill_viridis(trans = "log1p")+
  theme_classic()


# Extract the excretion data
  # This code summarized the excretion within each quad by averaging across the extrapolated values contained in each 1 m2 quadrat. 

#test the interspection function
temp2 <- df %>% group_by(plot, period, status, sample) %>%
  st_join(quads, join = st_intersects) %>%
  filter(plot == "CC1", period == "Closed season", sample == "A") %>%
  drop_na(transect, quad)

ggplot(temp2)+
  geom_sf(aes(geometry = geometry, fill =z)) # ok so it is extracting the data from under the quads!!!


extracted_data <- df %>% group_by(plot, period, status, sample) %>%
  st_join(quads, join = st_intersects) %>%
  drop_na(transect, quad) %>%
  ungroup() %>% 
  group_by(plot, period, status, transect, sample, quad) %>% 
  summarize(mean_nh4 = mean(z))

write_csv(extracted_data, file = "data/extracted_extretionestimates_byquad.csv")












