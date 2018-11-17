#this model looks at tail loss (1/0) as a function of the proportion Peaks of Otter
#salamanders at a site
#includes a random effect for site (propPO measured at this level) 
#and for pass/date as nested random effects
#so model is TLbin ~ propPO + (1|site) + (1|pass/date)

library(rstan)

dat<-podata_tl  #also podata_tl

#adds an integer index for survey date for stan
dat$date.int<-as.integer(as.factor(rbdata_tl$date2))

stan_code= '
data {
int<lower=1> N;                         //length of vector
int<lower=1> N_dates;                   //number of dates
int<lower=1> N_surveys;                 // number of surveys
int<lower=1,upper=N_surveys> survey[N];
int<lower=1,upper=N_dates> date[N];
int<lower=0,upper=1> y[N];  //tail loss
int<lower=0> n_site; 
int<lower=1,upper=n_site> site[N];
vector<lower=0,upper=1> [N] propPO;
}

parameters {
matrix[N_dates, N_surveys] a; // date and survey intercepts
vector[N_surveys] mu_survey;
vector[n_site] s;
vector[2] b;
real<lower=0> sigma_survey;
real<lower=0> sigma_date;
real<lower=0> sigma_site;

}
model {
vector[N] y_hat;
for (i in 1:N) 
y_hat[i] = b[1] + b[2]*propPO[i] + s[site[i]] + a[date[i], survey[i]]; 

//priors
sigma_survey ~ cauchy(0, 2.5);
sigma_date~cauchy(0, 2.5);
sigma_site~cauchy(0,2.5);
b~normal(0,3);
s ~ normal (0, sigma_site);

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
                            y = dat$tlbin,
                            n_site = 12,
                            N_dates = 30, 
                            N_surveys = 8,
                            propPO=dat$propPO,
                            site=dat$site.index,
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
apply(post,2, quantile, probs = c(0.025, 0.975))  #operations of the dataframe

stan_plot(stan_run, pars = c("sigma_survey", "sigma_date", "sigma_site"))
stan_plot(stan_run, pars = c("b[1]", "b[2]"))

print(stan_run, pars = c("b[1]", "b[2]"), digits_summary=3)
print(stan_run, pars = c("sigma_survey", "sigma_date", "sigma_site"), digits_summary=3)


#proportion of slope values from the posterior that are negative (decreased tail loss at sites with higher propPO)
mean(post$'b[2]'<0)



