
data {
  // data input
  int<lower=0> df_beta; // degrees of freedom of the observed difference in historical data
  int<lower=0> N; // total sample size for of the second stage (current data)
  vector<lower=0, upper=1>[N] group; // group indicator
  vector[N] y; // outcome (current data)
  
  // prior specificiations
  real diff_prior;
  real<lower=0> sigma_prior;
  real diff_skeptical;
  real<lower=0> sigma_skeptical;
  real<lower=0, upper=1> delta;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  target += log_mix(delta, student_t_lpdf(beta | df_beta, diff_prior, sigma_prior), student_t_lpdf(beta | 3, diff_skeptical, sigma_skeptical)); 
  target += normal_lpdf(y | alpha + beta * group, sigma);
}

