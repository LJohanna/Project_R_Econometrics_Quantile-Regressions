#################################
##   QUANTILE SELECTION MODEL  ##
#################################


###################
#  Preprocessing  #
###################

library(foreign)
library(quantreg)
if(!require(pacman)) install.packages("pacman")
pacman::p_load(copula)
library(copula)
library(gsl)
library(lpSolve)
library(lpSolveAPI)
library(stargazer)

# Load data
donnees <- read.dta("C:/Users/robin/Documents/ENSAE/3ème année/Semi and non-parametric econometrics/wagedata.dta")

# Remove missing values
donnees = donnees[which(donnees$wage != "NA") ,]
donnees = donnees[which(donnees$age != "NA") ,]
donnees = donnees[which(donnees$numhhkid != "NA") ,]
donnees = donnees[which(donnees$year != "NA") ,]
donnees = donnees[which(donnees$couple != "NA") ,]
donnees = donnees[which(donnees$ed16 != "NA") ,]
donnees = donnees[which(donnees$ed17 != "NA") ,]
donnees = donnees[which(donnees$ed18 != "NA") ,]
donnees = donnees[which(donnees$sex != "NA") ,]

# Remove negative wages
donnees = donnees[which(donnees$wage > 1) ,]


# We add noise to numbers because there is an error of singular design matrix coming from the duplicating 
# observations (for a single x value, multiple responses) 
donnees[,'wage'] = jitter(donnees[,'wage'])
donnees[,'age'] = jitter(donnees[,'age'])
donnees[,'numhhkid'] = jitter(donnees[,'numhhkid'])
donnees[,'year'] = jitter(donnees[,'year'])
donnees[,'couple'] = jitter(donnees[,'couple'])
donnees[,'ed16'] = jitter(donnees[,'ed16'])
donnees[,'ed17'] = jitter(donnees[,'ed17'])
donnees[,'ed18'] = jitter(donnees[,'ed18'])


# We separate males and females
donneeshommes = donnees[which(donnees$sex == 1) ,]
donneesfemmes = donnees[which(donnees$sex == 2) ,]

# Isolate workers and women
wageworkers = donneesfemmes[which(donneesfemmes$work == 1) ,] # Only select workers
N = nrow(donneesfemmes) # Number of data points
Nworkers = nrow(wageworkers) # Number of workers

#Quantile regression with tau=0.9
rq1 = rq(log(wage) ~ age + numhhkid + year + ed16 + ed17 + ed18, tau = 0.9,
         data=wageworkers)
summary(rq1)


#Quantile
print(mean(matrix(c(rep(1,Nworkers),wageworkers[,'age'],wageworkers[,'numhhkid'],wageworkers[,'year'],wageworkers[,'ed16'],wageworkers[,'ed17'],wageworkers[,'ed18']),nrow=Nworkers,ncol=7,byrow = FALSE)%*%rq1$coefficients))


#Quantile regression with tau=0.1
rq2 = rq(log(wage) ~ age + numhhkid + year  + ed16 + ed17 + ed18, tau = 0.1,
         data=wageworkers)
summary(rq2)

#Quantile
print(mean(matrix(c(rep(1,Nworkers),wageworkers[,'age'],wageworkers[,'numhhkid'],wageworkers[,'year'],wageworkers[,'ed16'],wageworkers[,'ed17'],wageworkers[,'ed18']),nrow=Nworkers,ncol=7,byrow = FALSE)%*%rq2$coefficients))


##############
# First step #
##############

# FIRST STEP : probit model with the following excluded variable (couple and ben_inc)

probit <- glm(work ~ age + numhhkid + year  + ed16 + ed17 + ed18 + couple + ben_inc,
              family=binomial(link="probit"),
              data=donneesfemmes)
summary (probit)


# We add a column to wageworkers in order to compute the constant in the quantile regressions
wageworkers['add'] = rep(1, Nworkers)

# We construct the vector of probabilities of working from the probit model for all individuals
propscoredf = subset(wageworkers, select=c('age','numhhkid','year','ed16','ed17','ed18','couple','ben_inc'))
propscore = predict(probit, propscoredf, type="response")
propscore

# Values for the grid search
copvalues = c(seq(from=0.1, to=9.9, by=0.2), seq(from=10, to=200, by=5), seq(from=-0.1, to=-9.9, by=-0.2))
length(copvalues)


#################
#  SECOND STEP  #
#################

vectau = c(runif(1), runif(1), runif(1), runif(1), runif(1), runif(1), runif(1), runif(1), runif(1), runif(1)) # Finite grid

ones = rep(1, Nworkers)
X = data.matrix(subset(wageworkers, select=c('add', 'age','numhhkid','year','ed16','ed17','ed18')))
Y = data.matrix(subset(wageworkers, select = c('wage')))
Y = log(Y) # Take log of the wage



# Compute objective function

objective_function = vector(length = length(copvalues)) # indexes allowing us to find the optimal rho
L = length(vectau)
C = length(copvalues)


for (i in 1:C) {
  
  copvalue = copvalues[i]
  objective = 0
  
  for (k in 1:L) {
 
    tau = vectau[k]
    
    # Frank copula
    frank_cop <- frankCopula(copvalue, dim=2)
    mx = matrix(c(rep(tau, Nworkers), propscore), 
                nrow=Nworkers, 
                ncol=2, byrow=FALSE)
    G_tau = matrix(pCopula(mx, frank_cop)/propscore)
    
    # Beta_tau
    fbeta = function(b){
      return(ones %*% (G_tau*max(Y - X%*%b, 0) + (1-G_tau)*max(X%*%b - Y, 0)))
    }
    # Optimization
    optim = optim(matrix(rep(0,7)), fbeta, method = "L-BFGS-B") 
    beta_tau <- optim$par
    
    # Objective function
    objective = objective + mean(propscore * (Y <= X %*% beta_tau) - G_tau)
    
  }
  
  objective_function[i] = objective^2 # Objective function is a norm
  print(i)
}


# Final objective function
objective_function

# Minimization of objective function
min_value = min(objective_function)
arg_min = which(objective_function==min_value)
rho = copvalues[arg_min]
rho

beta_tau


################
#  Third step  #
################

# Calculate the coefficients for other quantiles (in our case, tau = 0.1 and tau = 0.9)

# Tau = 0.1

# Frank copula
frank_cop <- frankCopula(rho, dim=2)
mx = matrix(c(rep(0.1, Nworkers), propscore), 
            nrow=Nworkers, 
            ncol=2, byrow=FALSE)
G_tau = matrix(pCopula(mx, frank_cop)/propscore)

# Beta_tau
fbeta = function(b){
  return(ones %*% (G_tau*pmax(Y - X%*%b, 0) + (1-G_tau)*pmax(X%*%b - Y, 0)))
}
# Optimization
optim_01 = optim(matrix(rep(0,7)), fbeta, method="L-BFGS-B") 
beta_tau_01 <- optim_01$par
beta_tau_01

# Tau = 0.9

# Frank copula
frank_cop <- frankCopula(rho, dim=2)
mx = matrix(c(rep(0.9, Nworkers), propscore), 
            nrow=Nworkers, 
            ncol=2, byrow=FALSE)
G_tau = matrix(pCopula(mx, frank_cop)/propscore)

# Beta_tau
fbeta = function(b){
  return(ones %*% (G_tau*pmax(Y - X%*%b, 0) + (1-G_tau)*pmax(X%*%b - Y, 0)))
}
# Optimization
optim_09 = optim(matrix(rep(0,7)), fbeta, method="L-BFGS-B") 
beta_tau_09 <- optim_09$par
beta_tau_09



############################################################################



# Isolate workers and men
wageworkers = donneeshommes[which(donneeshommes$work == 1) ,] # Only select workers
N = nrow(donneeshommes) # Number of data points
Nworkers = nrow(wageworkers) # Number of workers

#Quantile regression with tau=0.9
rq1 = rq(log(wage) ~ age + numhhkid + year + ed16 + ed17 + ed18, tau = 0.9,
         data=wageworkers)
summary(rq1)


#Quantile
print(mean(matrix(c(rep(1,Nworkers),wageworkers[,'age'],wageworkers[,'numhhkid'],wageworkers[,'year'],wageworkers[,'ed16'],wageworkers[,'ed17'],wageworkers[,'ed18']),nrow=Nworkers,ncol=7,byrow = FALSE)%*%rq1$coefficients))


#Quantile regression with tau=0.1
rq2 = rq(log(wage) ~ age + numhhkid + year  + ed16 + ed17 + ed18, tau = 0.1,
         data=wageworkers)
summary(rq2)

#Quantile
print(mean(matrix(c(rep(1,Nworkers),wageworkers[,'age'],wageworkers[,'numhhkid'],wageworkers[,'year'],wageworkers[,'ed16'],wageworkers[,'ed17'],wageworkers[,'ed18']),nrow=Nworkers,ncol=7,byrow = FALSE)%*%rq2$coefficients))


##############
# First step #
##############

# FIRST STEP : probit model with the following excluded variable (couple and ben_inc)

probit <- glm(work ~ age + numhhkid + year  + ed16 + ed17 + ed18 + couple + ben_inc,
              family=binomial(link="probit"),
              data=donneeshommes)
summary (probit)


# We add a column to wageworkers in order to compute the constant in the quantile regressions
wageworkers['add'] = rep(1, Nworkers)

# We construct the vector of probabilities of working from the probit model for all individuals
propscoredf = subset(wageworkers, select=c('age','numhhkid','year','ed16','ed17','ed18','couple','ben_inc'))
propscore = predict(probit, propscoredf, type="response")
propscore

# Values for the grid search
copvalues = c(seq(from=0.1, to=9.9, by=0.2), seq(from=10, to=200, by=5), seq(from=-0.1, to=-9.9, by=-0.2))
length(copvalues)


#################
#  SECOND STEP  #
#################

vectau = c(runif(1), runif(1), runif(1), runif(1), runif(1), runif(1)) # Finite grid

ones = rep(1, Nworkers)
X = data.matrix(subset(wageworkers, select=c('add', 'age','numhhkid','year','ed16','ed17','ed18')))
Y = data.matrix(subset(wageworkers, select = c('wage')))
Y = log(Y) # Take log of the wage



# Compute objective function

objective_function = vector(length = length(copvalues)) # indexes allowing us to find the optimal rho
L = length(vectau)
C = length(copvalues)


for (i in 1:C) {
  
  copvalue = copvalues[i]
  objective = 0
  
  for (k in 1:L) {
    
    tau = vectau[k]
    
    # Frank copula
    frank_cop <- frankCopula(copvalue, dim=2)
    mx = matrix(c(rep(tau, Nworkers), propscore), 
                nrow=Nworkers, 
                ncol=2, byrow=FALSE)
    G_tau = matrix(pCopula(mx, frank_cop)/propscore)
    
    # Beta_tau
    fbeta = function(b){
      return(ones %*% (G_tau*max(Y - X%*%b, 0) + (1-G_tau)*max(X%*%b - Y, 0)))
    }
    # Optimization
    optim = optim(matrix(rep(0,7)), fbeta, method = "L-BFGS-B") 
    beta_tau <- optim$par
    
    # Objective function
    objective = objective + mean(propscore * (Y <= X %*% beta_tau) - G_tau)
    
  }
  
  objective_function[i] = objective^2 # Objective function is a norm
  print(i)
}


# Final objective function
objective_function

# Minimization of objective function
min_value = min(objective_function)
arg_min = which(objective_function==min_value)
rho = copvalues[arg_min]
rho

beta_tau


################
#  Third step  #
################

# Calculate the coefficients for other quantiles (in our case, tau = 0.1 and tau = 0.9)


# Tau = 0.1

# Tau = 0.1

# Frank copula
frank_cop <- frankCopula(rho, dim=2)
mx = matrix(c(rep(0.1, Nworkers), propscore), 
            nrow=Nworkers, 
            ncol=2, byrow=FALSE)
G_tau = matrix(pCopula(mx, frank_cop)/propscore)

# Beta_tau
fbeta = function(b){
  return(ones %*% (G_tau*pmax(Y - X%*%b, 0) + (1-G_tau)*pmax(X%*%b - Y, 0)))
}
# Optimization
optim_01 = optim(matrix(rep(0,7)), fbeta, method="L-BFGS-B") 
beta_tau_01 <- optim_01$par
beta_tau_01



# Tau = 0.9

# Frank copula
frank_cop <- frankCopula(rho,dim=2)
mx = matrix(c(rep(0.1, Nworkers), propscore), 
            nrow=Nworkers, 
            ncol=2, byrow=FALSE)
G_tau = matrix(pCopula(mx, frank_cop)/propscore)

# Beta_tau
fbeta = function(b){
  return(ones %*% (G_tau*pmax(Y - X%*%b, 0) + (1-G_tau)*pmax(X%*%b - Y, 0)))
}
# Optimization
optim_09 = optim(matrix(rep(0,7)), fbeta, method="L-BFGS-B") 
beta_tau_09 <- optim_09$par
beta_tau_09



