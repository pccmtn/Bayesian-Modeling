library(rethinking)
data("WaffleDivorce")

################################################################################################
#############################              Load Dataset          ###############################
#############################             WaffleDivorce          ###############################
################################################################################################

d <- WaffleDivorce

################################################################################################
#############################              Spurious              ###############################
#############################             Associations           ###############################
################################################################################################

# standardize predictor

d$MedianAgeMarriage.s <- (d$MedianAgeMarriage - mean(d$MedianAgeMarriage)) / sd(d$MedianAgeMarriage)

# Fit bivariate model with age

m5.1 <- map(
  alist(
    Divorce ~ dnorm(mu, sigma),
    mu <- a + bA * MedianAgeMarriage.s,
    a ~ dnorm(10, 10),
    bA ~ dnorm(0,1),
    sigma ~ dunif(0,10)
  ), data = d
)

# compute percentile interval of mean

MAM.seq <- seq(from = -3, to = 3.5, length.out = 30)
mu <- link(m5.1, data.frame(MedianAgeMarriage.s = MAM.seq))
mu.PI <- apply(mu, 2, PI)

# plot

plot( Divorce ~ MedianAgeMarriage.s, data = d, col = rangi2)
abline(m5.1)
shade(mu.PI, MAM.seq)

# inspect precis output

precis(m5.1)

# each additional sd of delay in marriage (1.24 yrs) predicts a decrease of about one divorce per 1000
# adults (89% CI = -1.4, -0.7)

# Fit bivariate model with marriage rate

d$Marriage.s  <- (d$Marriage - mean(d$Marriage))/sd(d$Marriage)

m5.2 <- map(
  alist(
    Divorce ~ dnorm(mu, sigma),
    mu <- a + bM * Marriage.s,
    a ~ dnorm(10,10),
    bM ~ dnorm(0,1),
    sigma ~ dunif(0, 10)
  ), data = d
  )

MAM.seq <- seq(from = -4, to = 4.5, length.out = 30)
mu <- link(m5.2, data.frame(Marriage.s = MAM.seq))
mu.PI <- apply(mu, 2, PI)

plot(Divorce ~ Marriage.s, data = d, col = rangi2)
abline(m5.2)
shade(mu.PI, MAM.seq)

precis(m5.2)

################################################################################################
#############################           Multivariate             ###############################
#############################             Notation               ###############################
################################################################################################

m5.3 <- map(
  alist(
    Divorce ~ dnorm(mu, sigma),
    mu <- a + bR*Marriage.s + bA*MedianAgeMarriage.s,
    a ~ dnorm(10,10),
    bR ~ dnorm(0,1),
    bA ~dnorm(0,1),
    sigma ~ dunif(0,10)
  ),
  data = d)

precis(m5.3)
plot(precis(m5.3))

# once we know the median age at marriage for a state, there is little or no 
# additional predictive power in also knowing the rate of marriage in that state
# because the estimate for bR is so close to zero from both sides.

################################################################################################
#############################       Plotting Multivariate        ###############################
#############################             Posteriors             ###############################
################################################################################################

#### Predictor residual plots

# Predictor residual = the average prediction error when we use all of the other predictor variables to 
# model a predictor of interest.

# Leaves in the variation that is not expected by the model of the mu, as a function of the other
#predictors.

m5.4 <- map(
  alist(
    Marriage.s ~ dnorm(mu, sigma),
    mu <- a + b*MedianAgeMarriage.s,
    a ~dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~dunif(0,10)
  ),
  data = d)

### compute residuals by subtracting the observed rate from the predicted rate based on age at marriage

# compute expected value at MAP, for each state
mu <- coef(m5.4)['a'] + coef(m5.4)['b']*d$MedianAgeMarriage.s

# compute residual for each state

m.resid <- d$Marriage.s - mu

# plot

plot(Marriage.s ~ MedianAgeMarriage.s, d, col = rangi2)
abline(m5.4)

# loop over states
for (i in 1:length(m.resid)) {
  x <- d$MedianAgeMarriage.s[i] # x location of line segment
  y <- d$Marriage.s[i]
  #draw line
  lines(c(x,x), c(mu[i], y), lwd = 0.5, col= col.alpha("black", 0.7))
}

# these residuals are variation that is left over after taking out the purely linear relationship 
# between the two variables.

### How to use the residuals?

# Plot them against divorce rate, it will display the relationship between divorce and marriage rate
# after controlling for age at marriage. --> remaining association of each predictor with the outcome
# after already knowing the other predictors.

#### Posterior prediction plots

### Check how the model fits against the data

## Simulate prediction averaging over the posterior

# call link without specifying the new data, so it uses original data

mu <- link(m5.3)

# summarize samples across cases

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

# simulate new observations using original data

divorce.sim <- sim(m5.3, n = 1e4)
divorce.PI <- apply(divorce.sim, 2, PI)

# plot prediction against observed data

plot(mu.mean ~ d$Divorce, col = rangi2, ylim= range(mu.PI),
     xlab = "Observed Divorce", ylab = "Predicted Divorce")
abline(a = 0, b = 1, lty = 2)
for (i in 1:nrow(d)) {
  lines(rep(d$Divorce[i], 2), c(mu.PI[1, i], mu.PI[2,i]),
        col = rangi2)
}
identify(x = d$Divorce, y = mu.mean, labels = d$Loc, cex = 0.8)

# the model over predicts states with low divorce rates and under predicts high divorce rates

## Identify large model failures

# compute residuals

divorce.resid <- d$Divorce - mu.mean

# get ordering by divorce rate residuals

o <- order(divorce.resid)

# make the plot

dotchart(divorce.resid[o], labels = d$Loc[o], xlim = c(-6,5), cex = 0.6)
abline(v= 0, col = col.alpha("black", 0.2))
for ( i in 1:nrow(d)) {
  j <- o[i] # which state in order
  lines(d$Divorce[j] - c(mu.PI[1,j], mu.PI[2,j]), rep(i,2))
  points(d$Divorce[j] - c(divorce.PI[1,j], divorce.PI[2,j]), rep(i, 2),
         pch = 3, cex = 0.6, col= "grey")
}

## Final: plot residuals against new predictors to identify whether remaining variation is associated
# with another variable

waffles_per_capita <- d$WaffleHouses/d$Population
plot(divorce.resid ~ waffles_per_capita)
abline()

################################################################################################
#############################               Masked               ###############################
#############################             Relationships          ###############################
################################################################################################

data(milk)
d <- milk
str(d)

dcc <- d[complete.cases(d), ]

# To what extent energy content of milk is related to the percent of the brain mass?

#### First bivariate regression: kcal and neocortex

m5.5 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bn*neocortex.perc,
    a ~ dnorm(0,100),
    bn ~ dnorm(0,1),
    sigma ~ dunif(0,1)
  ), data = dcc)

precis(m5.5, digits = 3)

# a change from the smallest neocortex percent (55%) to the largest (76%), would result in an expected
# change in kcal of only:
coef(m5.5)["bn"]*(76-55)

# less than 0.1 kcalories

## Plot predicted mean and 89% interval for the mean

np.seq <- 0:100
pred.data <- data.frame(neocortex.perc = np.seq)

mu <- link(m5.5, data = pred.data, n = 1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(kcal.per.g ~ neocortex.perc, data = dcc, col = rangi2)
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1, ], lty = 2)
lines(np.seq, mu.PI[2,], lty = 2)

# weak positive relationship

# m5.6 is the same as m5.5 but using log body mass instead of neocortex
dcc$log.mass <- log(dcc$mass)


#### Multivariate regression

m5.7 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bn*neocortex.perc + bm*log.mass,
    a ~ dnorm(0,100),
    bn ~ dnorm(0,1),
    bm ~ dnorm(0,1),
    sigma ~ dunif(0,1)
  ), data = dcc)

precis(m5.7, digits = 3)
precis(m5.5, digits = 3) # see differences with bivariate model for bn


## Show how predicted kcal changes as a function of neocortex by keeping mass on hold on the mean

mean.log.mass <-(mean(dcc$log.mass))
np.seq <- 1:100

pred.data <- data.frame(
  neocortex.perc = np.seq,
  log.mass = mean.log.mass
)

mu <- link(m5.7, data = pred.data, n = 1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot(kcal.per.g ~ neocortex.perc, data = dcc, type = "n")
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1,], lty = 2)
lines(np.seq, mu.PI[2,], lty = 2)

################################################################################################
#############################                                    ###############################
#############################          Multicollinearity         ###############################
#############################                                    ###############################
################################################################################################

#### Multicollinear legs

# we try to predict height by using both left and right leg measurements

N <- 100
height <- rnorm(N, 10, 2)
leg_prop <- runif(N, 0.4, 0.5) # leg as proportion of height
leg_left <- leg_prop*height + rnorm(N, 0, 0.02) # sim left leg as prop of height + error
leg_right <- leg_prop*height + rnorm(N, 0, 0.02)

d <- data.frame(height, leg_left, leg_right)

m5.8 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left +br*leg_right,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2,10),
    br ~ dnorm(2,10),
    sigma ~ dunif(0,10)
  ), data = d
)
plot(precis(m5.8))

# if legs have almost identical lengths and height is so strongly associated with leg length,
# why does the prediction look so weird?

post <- extract.samples(m5.8)
plot(br ~ bl, post, col = col.alpha(rangi2, 0.1), pch = 16)

# the posterior distribution of the association of each leg with height is a narrow ridge of 
# negatively correlated values. The two parameters cannot be pulled apart and only their sum
# influences mu.

### simulating collinearity

if (FALSE) {
  d <- milk
  sim.coll <- function( r = 0.9) {
    d$x <- rnorm(nrow(d), mean = r*d$perc.fat,
                 sd = sqrt( (1-r^2)*var(d$perc.fat)))
    m <- lm(kcal.per.g ~ perc.fat + x, data = d)
    sqrt(diag(vcov(m)))[2]
  }
  rep.sim.coll <- function(r = 0.9, n = 100) {
    stddev <- replicate(n, sim.coll(r))
    mean(stddev)
  }
  
  r.seq <- seq(0, 0.99, by = 0.01)
  stddev <- sapply(r.seq, function(z) rep.sim.coll(r = z, n = 100))
  
  plot(stddev ~ r.seq, type = "l", col = rangi2, lwd = 2, xlab = "correlation")
}

################################################################################################
#############################                                    ###############################
#############################         Post-treatment bias        ###############################
#############################                                    ###############################
################################################################################################

N <- 100 #n of plants

h0 <- rnorm(N, 10, 2) # simulate initial heights

# assign treatments and simulate fungus and growth

treatment <- rep(0:1, each = N/2)
fungus <- rbinom(N, size = 1, prob = 0.5 - treatment*0.4)
h1 <- h0 + rnorm(N, 5 - 3*fungus)

# compose clean dataframe

d <- data.frame(h0 = h0, h1 = h1, treatment = treatment, fungus = fungus)

m5.13 <- map(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- a + bh*h0 + bt*treatment+bf*fungus,
    a ~ dnorm(0,100),
    c(bh, bt, bf) ~ dnorm(0,10),
    sigma ~ dunif(0,10)
  ), data = d
)

precis(m5.13)
plot(precis(m5.13))

# treatment has little to no effect, but we built the simulation in such a way that it matters to h1.
# what happened?

# fungus is a post-treatment variable because it is a consequence of treatment (treatment affects 
# growth by reducing fungus) --> exclude fungus

m5.14 <- map(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- a + bh*h0 + bt*treatment,
    a ~ dnorm(0,100),
    c(bh, bt) ~ dnorm(0,10),
    sigma ~ dunif(0,10)
  ), data = d
)

precis(m5.14)
plot(precis(m5.14)) # now treatment has a much bigger effect

################################################################################################
#############################                                    ###############################
#############################        Categorical Variables       ###############################
#############################                                    ###############################
################################################################################################

data("Howell1")
d <- Howell1
str(d)

m5.15 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bm*male,
    a ~ dnorm(178, 100),
    bm ~ dnorm(0, 10),
    sigma ~ dunif(0,50)
      ), data = d
)
precis(m5.15)
plot(precis(m5.15))

# a --> average height among females (assumes males = 0)
# bm --> average difference between males and females.

#### Many categories

d <- milk
unique(d$clade)
d$clade.NWM <- ifelse(d$clade == 'New World Monkey', 1, 0)
d$clade.OWM <- ifelse(d$clade == 'Old World Monkey', 1, 0)
d$clade.S <- ifelse(d$clade == 'Strepsirrhine', 1, 0)

# new model with dummies

m5.16 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + b.NWM*clade.NWM + b.OWM*clade.OWM+b.S*clade.S,
    a ~ dnorm(0.6, 10),
    c(b.NWM, b.OWM, b.S) ~ dnorm(0,1),
    sigma ~ dunif(0,10)
  ), data = d
)
precis(m5.16)
plot(precis(m5.16))

# a = average milk energy for Apes, parameters are differences between categories and the baseline Apes.

#### Other approach: Individual intercepts

d$clade_id <- coerce_index(d$clade)

m5.17 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a[clade_id],
    a[clade_id] ~ dnorm(0.6, 10),
    sigma ~dunif(0,10)
  ), data = d
)

precis(m5.17, depth = 2)
plot(precis(m5.17, depth = 2))
