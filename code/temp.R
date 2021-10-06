#--------------------------------------------------------------------------------
## Heat map on trial data for JP
# Bart DiFiore, August 2021
#--------------------------------------------------------------------------------

# The goal here is to estimate a statistical function describing the decline in nitrogen concentration as a function of lobster abundance and distance from den sites. Using this regression relationship, this code will then estimate the nitrogen concentration as a funciton of geographic space and lobster abundance. 


#---------------------------------------------------------------------------------
## Experimental regression relationship
#---------------------------------------------------------------------------------

df <- read.csv("Downloads/water_nutrients.csv") %>% janitor::clean_names()

ggplot(df, aes(x = distance, y = nh4_u_m))+
  geom_point(aes(size = count, color = site))

# Issue: What is your reponse (z) variable of interest? Extretion rates per indivdial, or water nitrogen concentrations? I ask because you have empirical estiamtes of nh4 at distance = 0 but you want to create a heat map of extretion rates ("EXCR_IND"). Is EXCR_IND in the same units? or different units. I'll build a heat map of nh4 based on the regression for now, and we can chat about how to convert between EXCR_IND. 

lm1 <- lm(log(nh4_u_m) ~ log(distance+0.001) + count, df)
summary(lm1)

newdat <- expand.grid(distance = seq(0,4,length.out = 100), count = mean(df$count))
newdat$nh4_u_m <- exp(predict(lm1, newdata = newdat))

ggplot(df, aes(x = distance, y = nh4_u_m))+
  geom_point(aes(size = count, color = site))+
  geom_line(data = newdat, aes(x = distance, y = nh4_u_m, group = count))
# scale_x_log10()+
# scale_y_log10()

# So this model feels adequate for now (R2 ~ 0.6). Certainly not perfect and we could refine this model, so that the initial decline isn't so steep using some other nonlinear functions. Essentially, it is a power law function such that N = 0.915 * distance ^ -0.15 + 0.05*count (I think? the second predictor might be throwing off my fast math). 

dist_decay <- function(distance, count, lm){
  log.y = coef(lm)[1] + coef(lm)[2]*log(distance+0.001) + coef(lm)[3]*count
  y = exp(log.y)
  y
} # Ignore this, I think the math is off somewhere here, and I used the predict() function rather than building my own.

#----------------------------------------------------------------------------------
## Extrapolate to one site
#----------------------------------------------------------------------------------

pl <- read.csv("Downloads/BR3_test.csv") %>% janitor::clean_names()

# for now I'm going to convert the plot data to counts at each coordinate, I'm doing this because the regression only uses distance and count to predict nitrogen. I don't think this is what you really want to do (becasue you want to model excr_ind, not nh4), but we can talk about how to convert between units and adapt the code to model excr_ind directly. 

ct <- pl %>% group_by(dist_x, dist_y) %>%
  dplyr::summarize(count = n())

ggplot(ct, aes(x = dist_x, y = dist_y))+
  geom_point(aes(size = count)) # This is a temporary plot displaying the count of lobsters are each unique location. 
# Issue: How did you determine if lobsters were at the same location? or at different locations? For instance some of the bugs at y = 0 in this plot, seem to be very close? How did you classify locaiton?


# Build a function to estimate nitrogen at every point in the grid. 

grid <- expand.grid(dist_x = seq(0, 20, length.out = 100), dist_y = seq(0,20, length.out = 100)) # construct a grid of points across which to estiamte nitrogen
plot(grid) # just confirm

points <- select(ct, dist_x, dist_y) # df of just the locations where lobsters were observed
count <- ct$count # vector of counts indexed by each point (in the same order as the points df above)

bd_fun <- function(points, grid){
  out <- t(fields::rdist(points, grid)) # estimate the euclidean distance from each focal point to each other point in the grid, transpose such that the matrix is organized where each row is a different grid point, and each column is a different focal point
  
  temp <- matrix(nrow = dim(out)[1], ncol = dim(out)[2]) # build a matrix to store the output of the for loop
  
  for(i in 1:dim(out)[2]){
    temp[,i] <- exp(predict(lm1, newdata = data.frame(distance = out[,i], count = count[i]))) # predict the contribution to the z variable from each focal point to each other point in the grid
  }
  
  z <- temp %*% matrix(data = 1, nrow = dim(out)[2], ncol = 1) # sum across rows in the matrix, such that it is the sum of the contributions to Z from each focal point where lobsters are located
  z # return Z
  
}

grid$z <- bd_fun(points = points, grid = grid) # run the function


ggplot(grid, aes(x = dist_x, y = dist_y))+ # plot it up! 
  geom_raster(aes(fill = z))+
  geom_point(data = points, aes(dist_x, dist_y, size = count), color = "red", pch = 21)+
  viridis::scale_fill_viridis()+
  labs(x = "", y = "", title = "Lobster sewage map")+
  theme_classic()

# Issue: So this plot (I think) is generally doing what you want it to. It is a deterministic estimate of how much nitrogen exists in every cell in the grid. Where the estimated nitrogen in a cell is the sum of the nitrogen contributions from every lobster den/location. The issue that I see currently is that I expected this plot to go to zero with distance from the site. As distance increases are are just starting to pick up background nitrogen production? Such that maybe we should scale the lower end of the distance decay function to zero? Also I'm still a bit stumped by the connection between EXCR_IND and nh4_u_m? Here I've predicted nh4_u_m based on the observed counts of lobsters but I have ignored the EXCR_IND variable in that data set. Chat more next week? 


# to visualize why the predictions are all > 20 here is a little toy example 

newdat <- expand.grid(distance = seq(0,100,length.out = 100), count = c(min(df$count), median(df$count), max(df$count)))
newdat$nh4_u_m <- exp(predict(lm1, newdata = newdat))

ggplot(newdat, aes(x = distance, y = nh4_u_m))+
  geom_line(aes(color = as.factor(count)))+
  coord_cartesian(ylim = c(0, 4)) # basically the log-log model that I built suggests that there is a lower asymptote (cause of the log-log scale), and that the number of lobsters at the focal site adjusts that lower asymptote up or down. One option would be to ignore lobster count and assume an addative effect. Then we could force the asymptote to approach 0, but I'd have to mull this one over a bit more.

# Maybe I'm overthinking this but one thought would be to estimate EXCR_IND for each of the locations where you have empirical data for nh4. Then we could use excr_ind as a predictor for observed nh4 (rather than count). That way we would be projecting nh4 as a function of excretion across space.
