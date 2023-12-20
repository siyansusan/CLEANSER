data {
  int<lower=1> N;
  array[N] int X;
  array[N] real L;
}

parameters {
  real<lower=0.000001> n_nbMean; // Noise negtive binomial mean rate parameter
  real<lower=0.1> n_nbDisp; //Noise negative binomial disp rate parameter
  real<lower=0.000001, upper=0.1> r;  // Mixture proportion
  real<lower=5> nbMean; // neg binomial mean parameter
  real<lower=0.1> nbDisp; // neg binomial dispersion parameter
}

model {
  // Priors:
  r ~ beta(1, 10);
  nbDisp ~ lognormal(log(3), log(2));
  n_nbMean ~ lognormal(log(1.075), log(1.1));
  n_nbDisp ~ lognormal(log(1.1), log(1.1));
  nbMean ~ lognormal(log(100), log(3));

  // Likelihoods:
  for(i in 1:N) {
    real sn_nbDisp = n_nbDisp;
    real sn_nbMean = n_nbMean * L[i];
    real snbDisp = nbDisp;
    real snbMean = nbMean * L[i];
    target += log_sum_exp(
                log(1-r) + neg_binomial_2_lpmf(X[i] | sn_nbMean, sn_nbDisp),
                log(r) + neg_binomial_2_lpmf(X[i] | snbMean, snbDisp)
              ) -
              log_diff_exp(
                log(1),
                log_sum_exp(
                  log(1-r) + sn_nbDisp * (log(sn_nbDisp) - log_sum_exp(log(sn_nbDisp), log(sn_nbMean))),
                  log(r) + snbDisp * (log(snbDisp) - log_sum_exp(log(snbMean), log(snbDisp)))
                )
              );
    }
}

generated quantities {
  array[N] real<lower=0,upper=1> PZi;
  array[N] real likeli_n_negbin;
  array[N] real likeli_negbin;

  for(i in 1:N) {
    real sn_nbDisp = n_nbDisp;
    real sn_nbMean = n_nbMean * L[i];
    real snbDisp = nbDisp;
    real snbMean = nbMean * L[i];
    real LP0 = log(1-r) + neg_binomial_2_lpmf(X[i] | sn_nbMean, sn_nbDisp);
    real LP1 = log(r) + neg_binomial_2_lpmf(X[i] | snbMean, snbDisp);
    PZi[i] = exp(LP1 - log_sum_exp(LP0, LP1));
    likeli_n_negbin[i] = LP0 - log(1-r);
    likeli_negbin[i] = LP1 - log(r);
  }
}

