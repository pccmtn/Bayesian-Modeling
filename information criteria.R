require(rethinking)

data("milk")

################################################################################################
#############################                 Using              ###############################
#############################         Information Criteria       ###############################
################################################################################################

d <- milk[complete.cases(milk), ]
d$neocortex <- d$neocortex.perc/100
dim(d)

a.start <- mean(d$kcal.per.g)
sigma.start <- log(sd(d$kcal.per.g))

# model with only intercept
m6.11 <- map(
  alist(
    kcal.per.g ~ dnorm(a, exp(log.sigma))
  ),
  data = d, start= list(a=a.start, log.sigma = sigma.start)
)

# model with only neocortex
m6.12 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, exp(log.sigma)),
    mu <- a +bn*neocortex
  ),
  data = d, start= list(a=a.start, bn = 0, log.sigma = sigma.start)
)

# model with only mass
m6.13 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, exp(log.sigma)),
    mu <- a +bm*log(mass)
  ),
  data = d, start= list(a=a.start, bm = 0, log.sigma = sigma.start)
)

# model with mass and neocortex
m6.14 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, exp(log.sigma)),
    mu <- a +bn*neocortex+bm*log(mass)
  ),
  data = d, start= list(a=a.start, bn = 0, bm = 0, log.sigma = sigma.start)
)

# all priors above are flat

### WAIC

WAIC(m6.14)

# first value == WAIC, smaller is better
# second value == lppd
# third value == pWAIC, for WAIC = -2(lppd - pWAIC)
# fourth value == se(WAIC) rough guidance of uncertainty in WAIC arising from sampling

### Compare models

(milk.models <- compare(m6.11, m6.12, m6.13, m6.14))
plot(milk.models, SE= T, dSE = T)

# filled dots: in sample deviance of each model
# open dots: WAIC
# lines SE of each WAIC
# triangles: SE of difference between each WAIC and top ranked WAIC

### Compare estimates

coeftab(m6.11, m6.12, m6.13, m6.14) #nobs == number of observations

plot(coeftab(m6.11, m6.12, m6.13, m6.14))

################################################################################################
#############################                                    ###############################
#############################            Model Averaging         ###############################
#############################                                    ###############################
################################################################################################

# Compute counterfactual predictions

# neocortex from 0.5 to 0.8

nc.seq <- seq(from = 0.5, to = 0.8, length.out = 30)

d.predict <- list(
  kcal.per.g = rep(0,30), #empty outcome
  neocortex = nc.seq,
  mass = rep(4.5, 30) # average mass
)
  
pred.m6.14 <- link(m6.14, data = d.predict)  
mu <- apply(pred.m6.14, 2, mean)  
  
mu.PI <- apply(pred.m6.14,2,PI)

# plot

plot(kcal.per.g ~ neocortex, d, col=rangi2)
lines(nc.seq,mu, lty=2)
lines(nc.seq,mu.PI[1,], lty=2)
lines(nc.seq,mu.PI[2,], lty=2)


### CREATE ENSEMBLE

milk.ensemble <- ensemble(m6.11, m6.12, m6.13, m6.14, data = d.predict)
mu <- apply(milk.ensemble$link, 2, mean)
mu.PI <- apply(milk.ensemble$link, 2, PI)

lines(nc.seq, mu)
shade(mu.PI, nc.seq)

# the solid regression line shows the average mu at each value on horizontal axis.
# the shaded area corresponds to the 89% percentile region.