###  R script for participants 

# Day 1 Session 2 - Hui-Walter models

# Two test - One population Hui-Walter models



#Load packages
library("tidyverse")
library("runjags")
library("rjags")
runjags.options(silent.jags=TRUE, silent.runjags=TRUE)
set.seed(2022-08-31)


### Two_test - One_population setting
#For demonstration purposes we start with the two_test - one_population setting (`hw_definition`) 


# The data are summarized in a two_x_two table (2^2 cells)
# or in a vector that contains all possible test results combinations

# Model Specification ('hw_definition')

hw_definition <- c("model{
  Cross_Classified_Data ~ dmulti(prob, N)
  
  # Test1+ Test2+
	prob[1] <- (prev * ((se[1])*(se[2]))) + ((1-prev) * ((1-sp[1])*(1-sp[2])))
  
  # Test1+ Test2-
	prob[2] <- (prev * ((se[1])*(1-se[2]))) + ((1-prev) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob[3] <- (prev * ((1-se[1])*(se[2]))) + ((1-prev) * ((sp[1])*(1-sp[2])))

  # Test1- Test2-
	prob[4] <- (prev * ((1-se[1])*(1-se[2]))) + ((1-prev) * ((sp[1])*(sp[2])))

  prev ~ dbeta(1, 1)
  se[1] ~ dbeta(1, 1)
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)
  sp[2] ~ dbeta(1, 1)

  #data# Cross_Classified_Data, N
  #monitor# prev, prob, se, sp
  #inits# prev, se, sp
}
")


# Specify - Load  model data
twoXtwo <- matrix(c(36, 4, 12, 48), ncol=2, nrow=2)
twoXtwo

Cross_Classified_Data <- as.numeric(twoXtwo)
N <- sum(Cross_Classified_Data)

# Initial values
prev <- list(chain1=0.05, chain2=0.95)
se <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))
sp <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))

# Run the model
results <- run.jags(basic_hw, n.chains=2)


# Remember to check convergence and effective sample size!


# Results interpretation and visualization
summary(results)


res <- summary(results)[,c(1:3,9,11)]
res[] <- round(res, 3)


plot(results)


pt <- plot(results)

print(pt[["prev.plot1"]])


print(pt[["se[1].plot1"]])
print(pt[["sp[1].plot1"]])


#### Label Switching
  
##  How to interpret a test with Se=0% and Sp=0%?
##  * The test is perfect - we are just holding it upside down...

## We can force Se+Sp >= 1.
## How? via the prior specification part in the model
    #  se[1] ~ dbeta(1, 1)
    #  sp[1] ~ dbeta(1, 1)T(1-se[1], )

        #Or:

    #  se[1] ~ dbeta(1, 1)T(1-sp[1], )
    # sp[1] ~ dbeta(1, 1)


# Try different prior distributions

# Weakly informative priors:
    # again via the prior specification part in the model

    # se[1] ~ dbeta(2, 1)
    # sp[1] ~ dbeta(2, 1)


# Stronger informative prior distributions

# A quick way to see the distribution of a prior:

# Beta(1,1)
curve(dbeta(x, 1, 1), from=0, to=1)
qbeta(c(0.025,0.975), shape1=1, shape2=1)

# Beta(2,1)
curve(dbeta(x, 2, 1), from=0, to=1)
qbeta(c(0.025,0.975), shape1=2, shape2=1)

## Choosing a prior

# What we want is e.g. Beta(20,1)

#But typically we have median and 95% confidence intervals from a paper, e.g.:

#"The median (95% CI) estimates of the sensitivity and specificity of the shiny new test 
# were 94% (92-96%) and 99% (97-100%) respectively"

# * How can we generate a Beta( , ) prior from this?

## The PriorGen package

# Median (95% CI) estimates of Se and Sp were 94% (92-96%) and 99% (97-100%)

library("PriorGen")
findbeta(themedian = 0.94, percentile=0.95, percentile.value = 0.92)


# Note: `themedian` could also be `themean`

# Plot beta distribution
curve(dbeta(x, shape1=429.95, shape2=27.76))


## Initial values

# Part of the problem before was also that we were specifying extreme initial values:

# se <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))
# sp <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))


# Let's change these to:
# se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
# sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))


## Exercise 1 
  
#  Run the `hw_definition` model under the following different scenarios and interpret the results in each case.

# 1. Change the priors for *Se* and *Sp* and try Beta(2,1).

# 2. Estimate the Beta parameters for the *Se* and *Sp* of the shiny new test that is described in the slides.

# 3. Run the model using the beta priors for *Se* and *Sp* from Step 2.

# 4. Try to run the model with different initial values.

# 5. Force Se + Sp > 1. Be careful to specify initial values that are within the restricted parameter space.



# Time for our lunch break



# Multi-population Hui-Walter models

## Hui-Walter models with multiple populations
  
#  - Basically an extension of the single-population model

# Assumptions

## 1. Different prevalence in different populations
  
# - In each population the data are summarized in a two_x_two table (2^2 cells) again
# - or each population has a vector that contains all possible test results combinations


## 2. The sensitivity and specificity must be consistent between populations

## 3. Account for conditional (in)dependence of the diagnostic tests (to be discussed tomorrow)

# Initial values
  
#  We have to be careful to make sure that the length of initial values for `prev` in each chain is equal to the number of populations

# For example with 5 populations we need:
  
    #prev <- list(chain1=c(0.1, 0.1, 0.1, 0.9, 0.9), chain2=c(0.9, 0.9, 0.9, 0.1, 0.1))

#The values you choose for different populations in the same chain can be the same - just make sure you pick different values for the same population between chains (i.e. *over-dispersed* initial values)

  ## Incorporating populations with known prevalence
  
#  Up to now prevalence has been a parameter, but it can also be (partially) observed:
  
# To fix the prevalence of population 1 we could do:
#Populations <- 5
#prev <- rep(NA, Populations)
#prev[1] <- 0
#prev

# But you also need to account for this in the initial values:
  # prev <- list(chain1=c(NA, 0.1, 0.1, 0.9, 0.9), chain2=c(NA, 0.9, 0.9, 0.1, 0.1))

## Data and initial value lists
  
  # There are actually multiple ways to specify data and initial values to runjags, including via the `data` and `inits` arguments

# Data
#data <- list(
#  Tally = Tally,
#  TotalTests = apply(Tally, 2, sum),
#  Populations = dim(Tally, 2),
#  prev = rep(NA, Populations),
# )
#data$prev[1] <- 0

# Initial values
#inits <- list(
#  chain1 = list(
#    prev = c(NA, 0.1, 0.1, 0.9, 0.9),
#    se = c(0.5, 0.99),
#    sp = c(0.5, 0.99)
#  ),
#  chain2 = list(
#    prev = c(NA, 0.9, 0.9, 0.1, 0.1),
#    se = c(0.99, 0.5),
#    sp = c(0.99, 0.5)
#  )
#)


# Run the modek
#results <- run.jags(..., data = data, inits = inits)


# Let's simulate a dataset to work on...

# Set a random seed so that the data are reproducible:
set.seed(2022-08-31)

sensitivity <- c(0.9, 0.6)
specificity <- c(0.95, 0.9)
N <- 1000

# Change the number of populations here:
Populations <- 5
# Change the variation in prevalence here:
(prevalence <- runif(Populations, min=0.1, max=0.9))

data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
  mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
  mutate(Test1 = rbinom(N, 1, sensitivity[1]*Status + (1-specificity[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status)))

(twoXtwoXpop <- with(data, table(Test1, Test2, Population)))
(Tally <- matrix(twoXtwoXpop, ncol=Populations))
(TotalTests <- apply(Tally, 2, sum))

# See data
Tally


## Exercise 2

# Start with 5 populations and analyse the data using the independent prevalence model.

# Now try ro run the same model with 2, and 3 populations, instead of 5.

    # - How does this affect the confidence intervals for the diagnostic test parameters?
  
# Now change the simulated prevalence so that it varies between 0.4-0.6 rather than 0.1-0.9.

    # - How does this affect the confidence intervals for the diagnostic test parameters?



# Hint - Help

## Model specification
  
multipop <- "
model{
  for(p in 1:Populations){
    Tally[1:4, p] ~ dmulti(prob[1:4, p], TotalTests[p])
  
    # Test1- Test2-
	 	  prob[1,p] <- (prev[p] * ((1-se[1])*(1-se[2]))) + ((1-prev[p]) * ((sp[1])*(sp[2])))
    # Test1+ Test2-
  	prob[2,p] <- (prev[p] * ((se[1])*(1-se[2]))) + ((1-prev[p]) * ((1-sp[1])*(sp[2])))
    # Test1- Test2+
  	prob[3,p] <- (prev[p] * ((1-se[1])*(se[2]))) + ((1-prev[p]) * ((sp[1])*(1-sp[2])))
  	 # Test1+ Test2+
  	 prob[4,p] <- (prev[p] * ((se[1])*(se[2]))) + ((1-prev[p]) * ((1-sp[1])*(1-sp[2])))

    prev[p] ~ dbeta(1, 1)
  }
  se[1] ~ dbeta(1, 1)T(1-sp[1], )
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)T(1-sp[2], )
  sp[2] ~ dbeta(1, 1)

  #data# Tally, TotalTests, Populations, se_prior, sp_prior
  #monitor# prev, se, sp
  #inits# prev, se, sp
  #module# lecuyer
}"


## Solution 2


# Set up initial values for 5 populations:
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
prev <- list(chain1=c(0.1, 0.1, 0.1, 0.9, 0.9), chain2=c(0.9, 0.9, 0.9, 0.1, 0.1))

# And run the model:
results_5p <- run.jags(multipop, n.chains=2)

# Remember to check convergence!
 plot(results_5p)
summary(results_5p)

# To change the number of populations and range of prevalence you just need to modify the simulation code, for example 3 populations with prevalence between 0.4-0.6 can be obtained using:

# Change the number of populations here:
Populations <- 3
# Change the variation in prevalence here:
(prevalence <- runif(Populations, min=0.4, max=0.6))

data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
  mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
  mutate(Test1 = rbinom(N, 1, sensitivity[1]*Status + (1-specificity[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status)))

(twoXtwoXpop <- with(data, table(Test1, Test2, Population)))
(Tally <- matrix(twoXtwoXpop, ncol=Populations))
(TotalTests <- apply(Tally, 2, sum))
# Adjust initial values for 3 populations:
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
prev <- list(chain1=c(0.1, 0.1, 0.9), chain2=c(0.9, 0.9, 0.1))


# And run the model:
results_3p <- run.jags(multipop, n.chains=2)
# Remember to check convergence!
plot(results_3p)
summary(results_3p)

# Note that when the effective sample size is not enough - you either need to run the model for longer in the first place, or extend it to get more samples:

# Extend the model:
results_3p <- extend.jags(results_3p, sample=50000)

# Remember to check convergence!
plot(results_3p)
summary(results_3p)

# As a general rule, the more populations you have, and the more the prevalence varies between them, the better.  However, this is conditional on having a consistent sensitivity and specificity between your populations!!!

