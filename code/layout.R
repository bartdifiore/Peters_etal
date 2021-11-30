library(tidyverse)
library(sp)
#----------------------------------------------------------------
## Extract excretion for analysis of benthic community structure
#----------------------------------------------------------------

coords <- read.csv("data/quadrat_coords.csv") %>% rename_all(tolower)

plot(dist_y ~ dist_x, coords)

utransects <- cbind(c(0, 5, 10, 15), c(0, 0, 0, 0), c(0, 5, 10, 15), c(20,20,20,20))
ID.trans <- c("A", "B", "C", "D")

a <- vector('list', length(2))

# loop through each centroid value and create a polygon
# this is where we match the ID to the new plot coordinates
for (i in 1:nrow(transects)) {  # for each for in object centroids
  a[[i]]<-sp::Lines(list(sp::Line(matrix(transects[i, ], ncol=2, byrow=TRUE))), ID.trans[i]) 
  # make it an Polygon object with the Plot_ID from object ID
}

# convert a to SpatialPolygon
transect.lines<-sp::SpatialLines(a)


# Now do the quadrats

ll <- cbind(coords$dist_x, coords$dist_y)
ul <- cbind(coords$dist_x, coords$dist_y + 1)
ur <- cbind(coords$dist_x + 1, coords$dist_y + 1)
lr <- cbind(coords$dist_x + 1, coords$dist_y )

ID <- paste(coords$transect, coords$quad, sep = "-")


square <- cbind(ll, ul, ur, lr, ll)


# First, initialize a list that will later be populated
# a, as a placeholder, since this is temporary
a <- vector('list', length(2))

# loop through each centroid value and create a polygon
# this is where we match the ID to the new plot coordinates
for (i in 1:nrow(coords)) {  # for each for in object centroids
  a[[i]]<-sp::Polygons(list(sp::Polygon(matrix(square[i, ], ncol=2, byrow=TRUE))), ID[i]) 
  # make it an Polygon object with the Plot_ID from object ID
}

# convert a to SpatialPolygon
polysB<-sp::SpatialPolygons(a)

sp::plot(polysB)

# build a line shape file for the transects




coordinates(coords) <- ~ dist_x + dist_y

png(filename = "figures/layout2.png", width = 1000, height = 750)
plot(coords)
plot(polysB, add = T)
text(coords, ID, cex = 0.5)
plot(transect.lines, add = T)
axis(side = c(1), at = c(-5, -1.5, 0, 1.5, 3.5, 5, 6.5, 8.5, 10, 11.5, 13.5, 15, 16.5, 20), labels = c(-5, -1.5, 0, 1.5, 3.5, 5, 6.5, 8.5, 10, 11.5, 13.5, 15, 16.5, 20))
axis(side = c(2))
dev.off()
