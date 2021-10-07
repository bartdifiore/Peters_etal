f(x) = c + (d-c)(\exp(-x/e))

exp.decay <- function(x, c, d, steepness){
  c + (d-c)*(exp(-x/steepness))
}

distance = 0:100

y <- exp.decay(x = distance, c = 0, d = 10, steepness = 10^1)
y2 <- exp.decay(x = distance, c = 0, d = 5, steepness = 10^1)

plot(y ~ distance)
points(y2 ~ distance, col = "red")




df <- read.csv("data/water_nutrients.csv",stringsAsFactors = F,na.strings=".")


mod1 <- nls(NH4_uM ~ c + (d-c)*(exp(-Distance/e)), data = df, start = list(c = 1, d = mean(df$NH4_uM[df$Distance == 0]), e = 2))
summary(mod1)

newdat <- data.frame(Distance = seq(0,4, length.out = 100), 
                     SFDM_g = mean(df$SFDM_g))
newdat$predicted <- predict(mod1, newdata = newdat)

ggplot(df, aes( x= Distance, y = NH4_uM))+
  geom_point(aes(size = SFDM_g))+
  geom_line(data = newdat, aes(x = Distance, y = predicted))



mll <- mle2(NH4_uM ~ dgamma(c + (d-c)*(exp(-Distance/e))),
             data=df,parameters = list(d~SFDM_g),
             start=list(c = 1, d = mean(df$NH4_uM[df$Distance == 0]), e = 2))
summary(mll)

newdat <- expand.grid(Distance = seq(0,4, length.out = 100), 
                     SFDM_g = c(min(df$SFDM_g),mean(df$SFDM_g), max(df$SFDM_g)))
newdat$predicted <- predict(mll, newdata = newdat)

ggplot(df, aes( x= Distance, y = NH4_uM))+
  geom_point(aes(size = SFDM_g))+
  geom_line(data = newdat, aes(x = Distance, y = predicted, group = SFDM_g))
