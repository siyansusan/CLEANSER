data {
  int<lower=1> N;
  array[N] int X;
  array[N] real L;
}
parameters {
  real<lower=0> lambda; // Poison rate parameter
  real<lower=0, upper=1.0> r;  // Mixture proportion
  real<lower=lambda> nbMean; // neg binomial mean parameter
  real<lower=0, upper=1.0> nbDisp; // neg binomial dispersion parameter
}
model {
   // Priors:
   r ~ beta(1, 10);
   lambda ~ lognormal(log(1.1), log(1.1));
   nbMean ~ lognormal(log(100), log(100));
   nbDisp ~ beta(1, 10);

   // Likelihoods:
   for(i in 1:N) {
      real slambda = lambda;
      real snbDisp = nbDisp;
      real snbMean = nbMean * L[i];
      target += log_sum_exp(log(1-r) + poisson_lpmf(X[i] | slambda),
                            log(r) + neg_binomial_2_lpmf(X[i] | snbMean, snbDisp)) -
                log_diff_exp(log(1),
                             log_sum_exp(log(1-r) + log(exp(-slambda)),
                                         log(r) + snbDisp * (log(snbDisp) - log_sum_exp(log(snbDisp),
                                                                                        log(snbMean)))));
      }
}
generated quantities {
   array[N]  real<lower=0,upper=1> PZi;
   array[N] real likeli_poisson;
   array[N] real likeli_negbin;
   real lambdasum=0;
   real nbsum=0;
   for(i in 1:N) {
      real slambda = lambda;
      real snbDisp = nbDisp;
      real snbMean = nbMean * L[i];
      real LP0 = log(1-r) + poisson_lpmf(X[i] | slambda);
      real LP1 = log(r) + neg_binomial_2_lpmf(X[i] | snbMean, snbDisp);
      PZi[i] = exp(LP1 - log_sum_exp(LP0, LP1));
      likeli_poisson[i] = LP0 - log(1-r);
      likeli_negbin[i] = LP1 - log(r);
   }
}

