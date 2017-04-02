library(rethinking)
data("Howell1")

################################################################################################
#############################         Load Dataset Howell        ###############################
#############################                                    ###############################
################################################################################################

d <- Howell1

d2 <- d[d$age >= 18, ]

dens(d2$height)

################################################################################################
#############################        Gaussian Model of           ###############################
#############################             Height                 ###############################
################################################################################################


curve(dnorm(x, 178,20), from = 100, to = 250)
curve(dunif(x, 0, 50), from = -10, to = 60)

sample_mu <- rnorm(1e4, 178, 10)
sample_sigma <- runif(1e4, 0, 60)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)


################################################################################################
#############################      Grid approximation of         ###############################
#############################  posterior distribution (6 steps)  ###############################
################################################################################################
  
# mu range of values
mu.list <- seq(140, 160, length.out = 200)

# sigma range of values
sigma.list <- seq(4,9, length.out = 200)

# matrix of all possible combinations of mu and sigma 
post <- expand.grid(mu = mu.list, sigma = sigma.list)

# Log likelihood of all possible combinations
post$LL <- sapply(1:nrow(post), function(i) sum( dnorm(
  d2$height,
  mean = post$mu[i], # compute individual LL and then sum() them all together 
  sd = post$sigma[i],
  log = T)))

# prior x LL
post$prod <- post$LL + dnorm(post$mu, 178, 20, T) + dunif(post$sigma, 0, 50, T)

# From log scale to probability scale by scaling by maximum log_product
post$prob <- exp(post$prod - max(post$prod))

# Graphical representations of results
image_xyz(post$mu, post$sigma, post$prob) # heatmap 
contour_xyz(post$mu, post$sigma, post$prob) # contour


################################################################################################
#############################         Sampling From the          ###############################
#############################       posterior distribution       ###############################
################################################################################################

sample.rows <- sample(1:nrow(post), size = 1e4, replace = T, prob = post$prob) # Randomly sample row number
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]
plot(sample.mu, sample.sigma, cex = 0.5, pch = 16, col = col.alpha(rangi2, 0.1)) # more plausible combinations of sigma and mu are where the plot is more dense

# marginal posteriors (averaging over the other parameter)
dens(sample.mu)
dens(sample.sigma)
HPDI(sample.mu)

################################################################################################
#############################      Quadratic approximation       ###############################
#############################       posterior distribution       ###############################
################################################################################################

# a-list of the definitions of the model 

flist <- alist(
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(156, 10),
  sigma ~ dunif(0,50)
)

# Fit model

m4.1 <- rethinking::map(flist, data = d2)

# Interpret results

precis(m4.1)

# The plausibility of each value of mu, after averaging over the plaudibilities of each value of sigma,
# is given by a gaussian with mean 154.6 and sd 0.4.
# The percentiles are interval boundaries of a 89% interval.

# visualize variance covariance matrix 
vcov(m4.1)

################################################################################################
#############################      Quadratic approximation       ###############################
#############################       posterior distribution       ###############################
#############################           using log_sigma          ###############################
################################################################################################

# Fit model

m4.1_logsigma <- rethinking::map(alist(
  height ~ dnorm(mu, log_sigma),
  mu ~ dnorm(156, 10),
  log_sigma ~ dnorm(2,10)
), data = d2)

# Interpret results

precis(m4.1_logsigma)

# The plausibility of each value of mu, after averaging over the plaudibilities of each value of sigma,
# is given by a gaussian with mean 154.6 and sd 0.4.
# The percentiles are interval boundaries of a 89% interval.

# extract samples and estimate sigma

post <- extract.samples(m4.1_logsigma)
sigma <- exp(post$log_sigma)


################################################################################################
#############################               Adding               ###############################
#############################            a Predictor             ###############################
################################################################################################

# We add weight and are interested in seeing how it covaries with height.

# First plot them so we have an idea about the relationship strength.

plot(d2$height ~ d2$weight)

# Fitting the model

m4.3 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(156, 100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ),
  data = d2)

# Interpreting the model fit

precis(m4.3) # gives the quadratic approximations of the parameters
# a = a person of 0 weight should have a height of 114
# b = a person 1 kg heavier is expected to be 0.90 cm taller
# sigma = width of the distribution of heights around the mean. 95% plausible heights lie within
# 10 cm (2sigma) of the mean of the height.


# Adding varcovar matrix

precis(m4.3, corr = TRUE)

# a and b are almost perfectly negatively correlated

# Solving correlation problem

# 1) Centering

d2$weight.c <- d2$weight - mean(d2$weight)

m4.4 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight.c,
    a ~ dnorm(156, 100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ),
  data = d2)

precis(m4.4, corr = TRUE)

# Now a means: the expected value of the outcome when the predictor is at its AVERAGE value.


################################################################################################
#############################               Plotting             ###############################
################################################################################################

plot(height ~ weight, data = d2)
abline(a = coef(m4.3)["a"], b = coef(m4.3)["b"])
# superimpose MAP values over height and weight


######################################################### Adding uncertainty around the mean

post <- extract.samples(m4.3)
post[1:5,] # the average of these lines is the MAP line

N <- 10 # 10 observations
dN <- d2[1:N,]
mN <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(156, 100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ),
  data = dN)

# extract 20 samples from the posterior 

post <- extract.samples(mN, n = 20)
plot(dN$weight, dN$height,
     xlim = range(dN$weight), ylim = range(dN$height),
     col= rangi2, xlab = "weight", ylab = "height")
mtext(concat("N = ", N))
for (i in 1:20) 
  abline(a = post$a[i], b = post$b[i], col = col.alpha("black", 0.3))

# vary N

N <- 100 # 10 observations
dN <- d2[1:N,]
mN <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(156, 100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ),
  data = dN)

# extract 20 samples from the posterior 

post <- extract.samples(mN, n = 20)
plot(dN$weight, dN$height,
     xlim = range(dN$weight), ylim = range(dN$height),
     col= rangi2, xlab = "weight", ylab = "height")
mtext(concat("N = ", N))
for (i in 1:20) 
  abline(a = post$a[i], b = post$b[i], col = col.alpha("black", 0.3))

plot1 <- ggplot(dN, aes(x = weight, y = height)) +
  geom_point() +
  xlim(range(dN$weight)) +
  ylim(range(dN$height))

## Adding uncertainty using contours

muat50 <- post$a + post$b*50 # vector of predicted means given weight = 50 kg.
#the variation across those means incorporates uncertainty in and correlation between both parameters.

HPDI(muat50, prob = 0.89)
# central 89% of the ways for the model to produce the data, place the average height between the upper
# and the lower bound of HPDI, assuming the weight is 50kg.

# we need to repeat this for every value of the weight --> Use Link function

mu <- link(m4.3)
str(mu)
#matrix of values of mu -- a distribution of mu for each individual in the original data

# we want to have a distribution for each individual weight value

weight.seq <- seq(from = 25, to = 70, by = 1) 

mu <- link(m4.3, data = data.frame(weight = weight.seq))
str(mu)

plot(height ~ weight, d2, type = "n")
for (i in 1:100)
  points(weight.seq, mu[i,], pch = 16, col = col.alpha(rangi2, 0.1))

# summarize the distribution of mu

mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.89)

# plot

plot(height ~ weight, data = d2, col = col.alpha(rangi2, alpha = 0.5))
lines(weight.seq, mu.mean)
shade(mu.HPDI, weight.seq)

################################################################################################
#############################           Link Function            ###############################
################################################################################################

post <- extract.samples(m4.3)
mu.link <- function(weight) post$a + post$b*weight # linear model of height
weight.seq <- seq(from = 25, to = 70, by = 1) # supply unique weights
mu <- sapply(weight.seq, mu.link) # s-apply link to each unique weight
mu.mean <- apply(mu, 2, mean) # compute mean of mu
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.89) # compute interval of mu


################################################################################################
#############################        Introducing sigma           ###############################
################################################################################################

sim.height <- sim(m4.3, data = list(weight = weight.seq)) #, n = 1e4)
str(sim.height)

height.PI <- apply(sim.height, 2, PI, prob = 0.89)

# plot

plot(height ~ weight, data = d2, col = col.alpha(rangi2, alpha = 0.5))

# draw map line
lines(weight.seq, mu.mean)

# draw HPDI region for line
shade(mu.HPDI, weight.seq)

# draw PI region for simulated weights
shade(height.PI, weight.seq)
