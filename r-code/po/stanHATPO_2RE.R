#model for proportion of hatchlings (compared to adults) for Peaks of Otter salamanders
#includes site type as a fixed effect (category) and pass/date nested random effects
#so model is stagebin ~ location + (1|pass/date)

library(rstan)

#prcoess data to add indicator variables for site type 
dat<-pohatch
summary(pohatch)
dat$dominant<-ifelse(dat$loc.int=='2', 1, 0)
dat$even<-ifelse(dat$loc.int=='3',1,0)
dat$rare<-ifelse(dat$loc.int=='4',1,0)


stan_code= '
data {
int<lower=1> N;                         //length of vector
int<lower=1> N_dates;                   //number of dates
int<lower=1> N_surveys;                 // number of surveys
int<lower=1,upper=N_surveys> survey[N];
int<lower=1,upper=N_dates> date[N];
int<lower=0,upper=1> y[N];              //stage class 0/1
vector<lower=0,upper=1>[N] dominant;
  vector<lower=0,upper=1>[N] even;
vector<lower=0,upper=1>[N] rare;
}

parameters {
matrix[N_dates, N_surveys] a; // date and survey intercepts
vector[N_surveys] mu_survey;
vector[4] b;
real<lower=0> sigma_survey;
real<lower=0> sigma_date;

}
model {
vector[N] y_hat;
for (i in 1:N) 
y_hat[i] = b[1] + b[2]*dominant[i] + b[3]*even[i] + b[4]*rare[i] + a[date[i], survey[i]]; 

sigma_survey ~ cauchy(0, 2.5);
sigma_date~cauchy(0, 2.5);
b~normal(0,3);

for (k in 1:N_surveys) { //loop through surveys 
  for (j in 1:N_dates) { //loop through dates 
    a[j,k] ~ normal (mu_survey[k], sigma_date);
    }
  mu_survey[k]~ normal(0, sigma_survey); 
}

y ~ bernoulli_logit(y_hat);
}
'

stan_run = stan(data = list(N = nrow(dat),
                            y = dat$stagebin,
                            N_dates = 31,  #number of survey dates
                            N_surveys = 8,  #number of passes
                            dominant=dat$dominant,
                            even=dat$even,
                            rare=dat$rare,
                            survey=dat$pass.int,
                            date=dat$date.int),  
                iter=4000,
                control = list(adapt_delta = 0.99, max_treedepth=15),
                model_code = stan_code)

#various summaries and plots
plot(stan_run)
summary(stan_run)  #summarizes parameter estimates
print(stan_run)   #includes measures of fit (rhat)
stan_dens(stan_run)  #nice density plots for posteriors

post<-as.data.frame(stan_run)   #extract parameter values into a dataframe
apply(post,2, quantile, probs = c(0.025, 0.975))  #operations on the dataframe

stan_plot(stan_run, pars = c("sigma_survey", "sigma_date"))
stan_plot(stan_run, pars = c("b[1]", "b[2]", "b[3]", "b[4]"))

print(stan_run, pars = c("b[1]", "b[2]", "b[3]", "b[4]"), digits_summary=3)
print(stan_run, pars = c("sigma_survey", "sigma_date"), digits_summary=3)

#proportion of posterior distribution below zero (i.e. fewer hatchlings compared to allopatric sites)
mean(post$'b[2]'<0)
mean(post$'b[3]'<0)
mean(post$'b[4]'<0)


