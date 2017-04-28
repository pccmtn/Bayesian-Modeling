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


