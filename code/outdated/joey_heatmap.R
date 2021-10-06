points <- data.frame(x = c(20,80), y = c(30,60)) # construct two locations of interest on the 100x100 grid

coords <- expand.grid(x = 1:100, y = 1:100) # build the grid

plot(y ~ x, coords)
points(y ~ x, points, col = "red", pch = 19)

dist <- seq(0,100, length.out = 100) # build a temporary distance vector to test the funciton
fun.z <- function(dist){ # build a sample distance decay funciton
  1*dist^-0.15
}
plot(fun.z(dist = dist)~dist) # plot the distance decay function over an arbitrary set of distances

out <- as.data.frame(t(fields::rdist(points, coords))) # estimate the euclidean distance from the two focal points of interest to every other point in the grid
names(out) <- c("dist.p1", "dist.p2") # rename the variables for teh distance from p1 and p2 (this is wonky and needs to be changes to sample n-points)


df <- cbind(coords, out) # put it all together in a data frame
df$z.p1 <- fun.z(dist = df$dist.p1) # estimate the response variable at each point in the grid as a function of distance from p1
df$z.p2 <- fun.z(dist = df$dist.p2) # estimate the response variable at each point in the grid as a function of distance from p2
df$z <- df$z.p1 + df$z.p2 # assuming that the observed concentration of z is additive, sum the contributions from each point source


ggplot(df, aes(x = x, y = y))+
  geom_raster(aes(fill = z))+
  geom_point(data = points, aes(x, y), color = "red", size = 2)+
  scale_fill_viridis()
