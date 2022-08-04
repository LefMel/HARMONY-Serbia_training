###  R script for participants 

# Day 2 Session 1 - Hui-Walter models accounting for conditional dependence

## Recap Hui-Walter Model assumptions 
# 1 The population is divided into two or more populations with different prevalences in which two or more tests are evaluated
# 2 Se and Sp are the same in all populations
# 3 The tests are *conditionally independent* given the disease status.

  ## Criticism of the Hui-Walter paradigm assumption (1)
  
# - The assumption (1) of distinct prevalences is necessary for the Hui-Walter model because otherwise, the data can be collapsed into a single 2x2 table with only three degrees of freedom for estimation.

# - The smaller the difference between disease prevalences, the larger are the posterior credible intervals, indicating a loss in precision.

# - The smallest difference in prevalence assessed by simulation was 10\%. In case of rare diseases, it might be difficult to find populations with prevalences higher than 10\%.

  ## Criticism of the Hui-Walter paradigm assumption (2)
  
#  - If assumption (2) is not satisfied, the accuracies would differ between two populations, this would add four additional parameters to be estimated, while there are only three additional degrees of freedom. 

# - *Se* and *Sp* are assumed to vary with external factors.

  ## Criticism of the Hui-Walter paradigm assumption (3)
  
#  - Assumption (3) demanding conditional independence was the first to be questioned by Vacek (1985).

# - Not accounting for potential conditional dependence may lead to misleading, biased estimates with a positive correlation leading to an over-estimation of the test accuracies and a negative of an under-estimation. 

  ### Condtional independencies
  
#  - Test are considered conditionally independent if the probability of getting a given test result on one test does not depend on the result from the other test, given the disease status of the individual.


# Example of a COVID-19 data set

## 2 pops - 3 tests

# Load packages
library("tidyverse")
library("runjags")
library("rjags")
runjags.options(silent.jags=TRUE, silent.runjags=TRUE)
set.seed(2022-09-01)


# Let's generate our dataset  

# The probability of infection with COVID in two populations:
prevalence <- c(0.01,0.05)
# The probability of shedding COVID in the nose conditional on infection:
nose_shedding <- 0.8
# The probability of shedding COVID in the throat conditional on infection:
throat_shedding <- 0.8
# The probability of detecting virus with the antigen test:
antigen_detection <- 0.75
# The probability of detecting virus with the PCR test:
pcr_detection <- 0.999
# The probability of random cross-reaction with the antigen test:
antigen_crossreact <- 0.05
# The probability of random cross-reaction with the PCR test:
pcr_crossreact <- 0.01

#  Simulating latent states:
N <- 20000
Populations <- length(prevalence)

covid_data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
  ## True infection status:
  mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
  ## Nose shedding status:
  mutate(Nose = Status * rbinom(N, 1, nose_shedding)) %>%
  ## Throat shedding status:
  mutate(Throat = Status * rbinom(N, 1, throat_shedding))


#  Simulating test results:
covid_data <- covid_data %>%
  ## The nose swab antigen test may be false or true positive:
  mutate(NoseAG = case_when(
    Nose == 1 ~ rbinom(N, 1, antigen_detection),
    Nose == 0 ~ rbinom(N, 1, antigen_crossreact)
  )) %>%
  ## The throat swab antigen test may be false or true positive:
  mutate(ThroatAG = case_when(
    Throat == 1 ~ rbinom(N, 1, antigen_detection),
    Throat == 0 ~ rbinom(N, 1, antigen_crossreact)
  )) %>%
  ## The PCR test may be false or true positive:
  mutate(ThroatPCR = case_when(
    Throat == 1 ~ rbinom(N, 1, pcr_detection),
    Throat == 0 ~ rbinom(N, 1, pcr_crossreact)
  ))


#  The overall sensitivity of the tests can be calculated as follows:
covid_sensitivity <- c(
  # Nose antigen:
  nose_shedding*antigen_detection + (1-nose_shedding)*antigen_crossreact,
  # Throat antigen:
  throat_shedding*antigen_detection + (1-throat_shedding)*antigen_crossreact,
  # Throat PCR:
  throat_shedding*pcr_detection + (1-throat_shedding)*pcr_crossreact
)
covid_sensitivity


# The overall specificity of the tests is more straightforward:
covid_specificity <- c(
  # Nose antigen:
  1 - antigen_crossreact,
  # Throat antigen:
  1 - antigen_crossreact,
  # Throat PCR:
  1 - pcr_crossreact
)
covid_specificity


## Model specification - snip
# prob[1,p] <-  prev[p] * ((1-se[1])*(1-se[2])*(1-se[3]) 
#                         +covse12 +covse13 +covse23) +
#  (1-prev[p]) * (sp[1]*sp[2]*sp[3] 
#                 +covsp12 +covsp13 +covsp23)

#prob[2,p] <- prev[p] * (se[1]*(1-se[2])*(1-se[3]) 
#                        -covse12 -covse13 +covse23) +
#  (1-prev[p]) * ((1-sp[1])*sp[2]*sp[3] 
#                 -covsp12 -covsp13 +covsp23)

## snip ##

# Covariance in sensitivity between tests 1 and 2:
#covse12 ~ dunif( (se[1]-1)*(1-se[2]) , 
#                 min(se[1],se[2]) - se[1]*se[2] )
# Covariance in specificity between tests 1 and 2:
#covsp12 ~ dunif( (sp[1]-1)*(1-sp[2]) , 
#                 min(sp[1],sp[2]) - sp[1]*sp[2] )

# It is quite easy to get the covariance terms slightly wrong!


  ## Template Hui-Walter function in runjags
  
#  The model code and data format for an arbitrary number of populations (and tests) can be determined automatically using the template_huiwalter function from the runjas package:
template_huiwalter(
  covid_data %>% select(Population, NoseAG, ThroatAG, ThroatPCR), 
  outfile = 'covidmodel.txt', covariance=TRUE)

# This generates self-contained model/data/initial values etc

# And can be run directly from R:
results <- run.jags('covidmodel.txt')
results

res <- summary(results)[,c(1:3,9,11)]
res[] <- round(res, 3)

---
  
  ## Template Hui-Walter
  
#  - Modifying priors must still be done directly in the model file
# - The model needs to be re-generated if the data changes
# * But remember that your modified priors will be reset

# - There must be a single column for the population (as a factor), and all of the other columns (either factor, logical or numeric) are interpreted as being test results

# - Covariance terms are all deactivated by default

  ## Activating covariance terms
  
#  Find the lines for the covariances that we want to activate (i.e. the two Throat tests):
  
#  You will also need to uncomment out the relevant initial values for BOTH chains

# There is also an argument in the template_huiwalter() function called covariance. When covariance=TRUE all covariance terms are automatically activated.
  
## Exercise 1
  
#  Use the template_huiwalter function() to look at the simple 2-test 5-population example from yesterday's session.  

# There are 3 steps for this exercise, specified below.

# Step 1 - simulate the data |  use this data simulation code:

# Set a random seed so that the data are reproducible:
set.seed(2022-09-01)
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
  mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status))) %>%
  select(-Status)

(twoXtwoXpop <- with(data, table(Test1, Test2, Population)))
(Tally <- matrix(twoXtwoXpop, ncol=Populations))
(TotalTests <- apply(Tally, 2, sum))

# Step 2 - run the model without accounting for conditional dependence between tests
template_huiwalter(data, outfile="template_2test.txt")


# Look at the model code and familiarise yourself with how the model is set out (there are some small differences, but the overall code is equivalent).Run the model.

# Step 3 - run the model accounting for conditional dependence between tests

# Now activate the correlation terms between tests 1 and 2.  Is anything different about the results?

# Step 4 (optional) - try with 2 populations - What's happens in this case? Is the model identifiable?

## Solution 1

# There is no particular solution to the first part of this exercise, but please ask if you have any questions about the model code that template_huiwalter generates.  Remember that re-running the template_huiwalter function will over-write your existing model including any changes you made, so be careful!

# We can run the model as follows:

results_nocov <- run.jags("template_2test.txt")
results_nocov

# A shortcut for activating the covariance terms is to re-run template_huiwalter as follows:
template_huiwalter(data, outfile="template_2test_cov.txt", covariance=TRUE)
results_cov <- run.jags("template_2test_cov.txt")
results_cov

