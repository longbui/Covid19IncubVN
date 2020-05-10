#authoir longbui189@gmail.com
#longbui189@gmail.com
#based on code of Backer et al. (https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.5.2000062)

#load library
library(tidyverse)
library(rstan)
library(loo)
#library(lubridate)
#stan config
options(auto_write = TRUE)
options(mc.cores = 4)
set.seed(12345)

#create input data
data.stan <- incub %>% 
  mutate(tSymptomOnset = SL,
         tStartExposure = EL,
         tEndExposure = ER) %>% dplyr::select(id, tSymptomOnset, tStartExposure,tEndExposure)

input.data <- list(
  N = nrow(data.stan),
  tStartExposure = data.stan$tStartExposure,
  tEndExposure = data.stan$tEndExposure,
  tSymptomOnset = data.stan$tSymptomOnset)


# compile model Weibull

model <- stan(data = input.data, 
              chains = 0, 
              iter = 0,
              model_code = "
data{
  int <lower = 1> N;
  vector[N] tStartExposure;
  vector[N] tEndExposure;
  vector[N] tSymptomOnset;
}

parameters{
  real<lower = 0> alphaInc; 	// Shape parameter of weibull distributed incubation period
  real<lower = 0> sigmaInc; 	// Scale parameter of weibull distributed incubation period
  vector<lower = 0, upper = 1>[N] uE;	// Uniform value for sampling between start and end exposure
}

transformed parameters{
  vector[N] tE; 	// moment of infection
  tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
}

model{
  // Contribution to likelihood of incubation period
  target += weibull_lpdf(tSymptomOnset -  tE  | alphaInc, sigmaInc);
}

generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = weibull_lpdf(tSymptomOnset[i] -  tE[i]  | alphaInc, sigmaInc);
  }
}
"
)
#fit stan, add control params for ensuring convergence

stanfit <- stan(fit = model, data = input.data, 
                init = "random",
                warmup = 4000,
                iter = 14000, 
                chains = 4,
                control = list(adapt_delta = 0.99))

print(stanfit) # check results and convergence
LL = extract_log_lik(stanfit, parameter_name = 'log_lik') # modelfit
loo1 <- loo(LL)
loo1

# results
alpha <- rstan::extract(stanfit)$alphaInc
sigma <- rstan::extract(stanfit)$sigmaInc

quantile(sigma*gamma(1+1/alpha), probs = c(0.025,0.5,0.975)) # posterior median and 95%CI of mean
quantile(sqrt(sigma^2*(gamma(1+2/alpha)-(gamma(1+1/alpha))^2)), probs = c(0.025,0.5,0.975)) # posterior median and 95%CI of sd


percentiles <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qweibull(p = p, shape = alpha, scale = sigma), probs = c(0.025, 0.5, 0.975)))
colnames(percentiles) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
percentiles # posterior median and 95%CI of percentiles

# Gamma distribution
model.gamma <- stan(data = input.data, 
              chains = 0, 
              iter = 0,
              model_code = "
data{
  int <lower = 1> N;
  vector[N] tStartExposure;
  vector[N] tEndExposure;
  vector[N] tSymptomOnset;
}

parameters{
  real<lower = 0> alphaG; 	// Shape parameter of weibull distributed incubation period
  real<lower = 0> betaG; 	// Rate parameter of weibull distributed incubation period
  vector<lower = 0, upper = 1>[N] uE;	// Uniform value for sampling between start and end exposure
}

transformed parameters{
  vector[N] tE; 	// infection moment
  tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
}

model{
  // Contribution to likelihood of incubation period
  target += gamma_lpdf(tSymptomOnset -  tE  | alphaG, betaG);
}

generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = gamma_lpdf(tSymptomOnset[i] -  tE[i]  | alphaG, betaG);
  }
}
"
)

#Gamma distribution fitting
stanfitG <- stan(fit = model.gamma, data = input.data, 
                init = "random",
                warmup = 4000,
                iter = 14000, 
                chains = 4, control = list(adapt_delta = 0.99))
print(stanfitG) # check results and convergence


LLG = extract_log_lik(stanfitG, parameter_name = 'log_lik') # Pareto Loo PSIS
loo2 <- loo(LLG)
loo2
alphaG <- rstan::extract(stanfitG)$alphaG #extract parameters
betaG <- rstan::extract(stanfitG)$betaG

quantile((alphaG / betaG), probs = c(0.025,0.5,0.975))# posterior median and 95%CI of mean
quantile(sqrt(alphaG / betaG^2), probs = c(0.025,0.5,0.975)) # posterior median and 95%CI of sd

percentilesG <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qgamma(p = p, shape = alphaG, scale = 1/betaG), probs = c(0.025, 0.5, 0.975)))
colnames(percentilesG) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
percentilesG  # posterior median and 95%CI of percentiles

#Log normal distribution fitting
model.log <- stan(data = input.data, 
                    chains = 0, 
                    iter = 0,
                    model_code = "
data{
  int <lower = 1> N;
  vector[N] tStartExposure;
  vector[N] tEndExposure;
  vector[N] tSymptomOnset;
}

parameters{
  real<lower = 0> mu; 	// mu parameter of lognormal distributed incubation period
  real<lower = 0> sigmaL; 	// sigma parameter of lognormal distributed incubation period
  vector<lower = 0, upper = 1>[N] uE;	// Uniform value for sampling between start and end exposure
}

transformed parameters{
  vector[N] tE; 	// infection moment
  tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
}

model{
  // Contribution to likelihood of incubation period
  target += lognormal_lpdf(tSymptomOnset -  tE  | mu, sigmaL);
}

generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = lognormal_lpdf(tSymptomOnset[i] -  tE[i]  | mu, sigmaL);
  }
}
"
)

stanfitL <- stan(fit = model.log, data = input.data, 
                init = "random",
                warmup = 4000,
                iter = 14000, 
                chains = 4, control = list(adapt_delta = 0.99))

LLL = extract_log_lik(stanfitL, parameter_name = 'log_lik')
loo3 <- loo(LLL)
loo3

mu <- rstan::extract(stanfitL)$mu
sigmaL <- rstan::extract(stanfitL)$sigmaL


quantile(exp(mu + (sigmaL^2)/2), probs = c(0.025,0.5,0.975)) # posterior median and 95%CI of mean
quantile(sqrt((exp(sigmaL^2)-1) * exp(2*mu + sigmaL^2)), probs = c(0.025,0.5,0.975))
percentilesL <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qlnorm(p = p, meanlog = mu, sdlog = sigmaL), probs = c(0.025, 0.5, 0.975)))
colnames(percentilesL) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
percentilesL

#loo comparing
sumW <- summary(stanfit)
sumG <- summary(stanfitG)
sumL <- summary(stanfitL)

loo_compare(loo1, loo2, loo3)


sumW$summary
tE <- summary(stanfit, pars = c("tE"))$summary
tE <- as.data.frame(tE)   

#summary(stanfit, pars = c("alphaInc", "sigmaInc"))$summary




              