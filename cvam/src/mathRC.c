/* ------------------------------------------------------------------ */
/* C functions callable from Fortran that access R's API ------------ */
/* ------------------------------------------------------------------ */
#include <R.h>
#include <Rmath.h>

/* ------------------------------------------------------------------ */
/* random number generator state ------------------------------------ */

void F77_SUB(getrandgenstaterc)(void) { GetRNGstate(); }

void F77_SUB(putrandgenstaterc)(void) { PutRNGstate(); }

/* ------------------------------------------------------------------ */
/* uniform distribution */

/* density function */
void F77_SUB(dunifrc)(double *x, double *a, double *b, 
		      int *logd, double *ans){
  double x_local = x[0];
  double a_local = a[0];
  double b_local = b[0];
  int logd_local = logd[0];
  ans[0] = dunif(x_local, a_local, b_local, logd_local);
}

/* distribution function */
void F77_SUB(punifrc)(double *q, double *a, double *b, 
		      int *lowertail, int *logp,  double *ans){
  double q_local = q[0];
  double a_local = a[0];
  double b_local = b[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = punif(q_local, a_local, b_local,
		 lowertail_local, logp_local);
}

/* quantile function */
void F77_SUB(qunifrc)(double *p, double *a, double *b, 
		      int *lowertail, int *logp,  double *ans){
  double p_local = p[0];
  double a_local = a[0];
  double b_local = b[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = qunif(p_local, a_local, b_local,
		 lowertail_local, logp_local);
}

/* random generator */
void F77_SUB(runifrc)(double *a, double *b, double *ans){
  double a_local = a[0];
  double b_local = b[0];
  ans[0] = runif(a_local, b_local);
}

/* ------------------------------------------------------------------ */
/* normal distribution */

/* density function */
void F77_SUB(dnormrc)(double *x, double *mu, double *sigma, 
		      int *logd, double *ans){
  double x_local = x[0];
  double mu_local = mu[0];
  double sigma_local = sigma[0];
  int logd_local = logd[0];
  ans[0] = dnorm(x_local, mu_local, sigma_local, logd_local);
}

/* distribution function */
void F77_SUB(pnormrc)(double *q, double *mu, double *sigma, 
		      int *lowertail, int *logp,  double *ans){
  double q_local = q[0];
  double mu_local = mu[0];
  double sigma_local = sigma[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = pnorm(q_local, mu_local, sigma_local,
		 lowertail_local, logp_local);
}

/* quantile function */
void F77_SUB(qnormrc)(double *p, double *mu, double *sigma, 
		      int *lowertail, int *logp,  double *ans){
  double p_local = p[0];
  double mu_local = mu[0];
  double sigma_local = sigma[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = qnorm(p_local, mu_local, sigma_local,
		 lowertail_local, logp_local);
}

/* random generator */
void F77_SUB(rnormrc)(double *mu, double *sigma, double *ans){
  double mu_local = mu[0];
  double sigma_local = sigma[0];
  ans[0] = rnorm(mu_local, sigma_local);
}

/* ------------------------------------------------------------------ */
/* chisquare distribution */

/* density function */
void F77_SUB(dchisqrc)(double *x, double *df, int *logd, double *ans){
  double x_local = x[0];
  double df_local = df[0];
  int logd_local = logd[0];
  ans[0] = dchisq(x_local, df_local, logd_local);
}

/* distribution function */
void F77_SUB(pchisqrc)(double *q, double *df,
		      int *lowertail, int *logp,  double *ans){
  double q_local = q[0];
  double df_local = df[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = pchisq(q_local, df_local, lowertail_local, logp_local);
}

/* quantile function */
void F77_SUB(qchisqrc)(double *p, double *df,
		      int *lowertail, int *logp,  double *ans){
  double p_local = p[0];
  double df_local = df[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = qchisq(p_local, df_local, lowertail_local, logp_local);
}

/* random generator */
void F77_SUB(rchisqrc)(double *df, double *ans){
  double df_local = df[0];
  ans[0] = rchisq(df_local);
}

/* ------------------------------------------------------------------ */
/* Student's t distribution */

/* density function */
void F77_SUB(dtrc)(double *x, double *df, int *logd, double *ans){
  double x_local = x[0];
  double df_local = df[0];
  int logd_local = logd[0];
  ans[0] = dt(x_local, df_local, logd_local);
}

/* distribution function */
void F77_SUB(ptrc)(double *q, double *df,
		      int *lowertail, int *logp,  double *ans){
  double q_local = q[0];
  double df_local = df[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = pt(q_local, df_local, lowertail_local, logp_local);
}

/* quantile function */
void F77_SUB(qtrc)(double *p, double *df,
		      int *lowertail, int *logp,  double *ans){
  double p_local = p[0];
  double df_local = df[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = qt(p_local, df_local, lowertail_local, logp_local);
}

/* random generator */
void F77_SUB(rtrc)(double *df, double *ans){
  double df_local = df[0];
  ans[0] = rt(df_local);
}

/* ------------------------------------------------------------------ */
/* beta distribution */

/* density function */
void F77_SUB(dbetarc)(double *x, double *shape1, double *shape2, 
		      int *logd, double *ans){
  double x_local = x[0];
  double shape1_local = shape1[0];
  double shape2_local = shape2[0];
  int logd_local = logd[0];
  ans[0] = dbeta(x_local, shape1_local, shape2_local, logd_local);
}

/* distribution function */
void F77_SUB(pbetarc)(double *q, double *shape1, double *shape2, 
		      int *lowertail, int *logp,  double *ans){
  double q_local = q[0];
  double shape1_local = shape1[0];
  double shape2_local = shape2[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = pbeta(q_local, shape1_local, shape2_local,
		 lowertail_local, logp_local);
}

/* quantile function */
void F77_SUB(qbetarc)(double *p, double *shape1, double *shape2, 
		      int *lowertail, int *logp,  double *ans){
  double p_local = p[0];
  double shape1_local = shape1[0];
  double shape2_local = shape2[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = qbeta(p_local, shape1_local, shape2_local,
		 lowertail_local, logp_local);
}

/* random generator */
void F77_SUB(rbetarc)(double *shape1, double *shape2, double *ans){
  double shape1_local = shape1[0];
  double shape2_local = shape2[0];
  ans[0] = rbeta(shape1_local, shape2_local);
}

/* ------------------------------------------------------------------ */
/* gamma distribution */

/* density function */
void F77_SUB(dgammarc)(double *x, double *shape, double *scale, 
		      int *logd, double *ans){
  double x_local = x[0];
  double shape_local = shape[0];
  double scale_local = scale[0];
  int logd_local = logd[0];
  ans[0] = dgamma(x_local, shape_local, scale_local, logd_local);
}

/* distribution function */
void F77_SUB(pgammarc)(double *q, double *shape, double *scale, 
		      int *lowertail, int *logp,  double *ans){
  double q_local = q[0];
  double shape_local = shape[0];
  double scale_local = scale[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = pgamma(q_local, shape_local, scale_local,
		 lowertail_local, logp_local);
}

/* quantile function */
void F77_SUB(qgammarc)(double *p, double *shape, double *scale, 
		      int *lowertail, int *logp,  double *ans){
  double p_local = p[0];
  double shape_local = shape[0];
  double scale_local = scale[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = qgamma(p_local, shape_local, scale_local,
		 lowertail_local, logp_local);
}

/* random generator */
void F77_SUB(rgammarc)(double *shape, double *scale, double *ans){
  double shape_local = shape[0];
  double scale_local = scale[0];
  ans[0] = rgamma(shape_local, scale_local);
}

/* ------------------------------------------------------------------ */
/* binomial distribution */

/* mass function */
void F77_SUB(dbinomrc)(double *x, double *n, double *prob, 
		      int *logd, double *ans){
  double x_local = x[0];
  double n_local = n[0];
  double prob_local = prob[0];
  int logd_local = logd[0];
  ans[0] = dbinom(x_local, n_local, prob_local, logd_local);
}

/* mass function continuous version */
void F77_SUB(dbinomrawrc)(double *x, double *n, double *probp, 
			  double *probq, int *logd, double *ans){
  double x_local = x[0];
  double n_local = n[0];
  double probp_local = probp[0];
  double probq_local = probq[0];
  int logd_local = logd[0];
  ans[0] = dbinom_raw(x_local, n_local, probp_local, probq_local, logd_local);
}

/* distribution function */
void F77_SUB(pbinomrc)(double *q, double *n, double *prob, 
		      int *lowertail, int *logp,  double *ans){
  double q_local = q[0];
  double n_local = n[0];
  double prob_local = prob[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = pbinom(q_local, n_local, prob_local,
		 lowertail_local, logp_local);
}

/* quantile function */
void F77_SUB(qbinomrc)(double *p, double *n, double *prob, 
		      int *lowertail, int *logp,  double *ans){
  double p_local = p[0];
  double n_local = n[0];
  double prob_local = prob[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = qbinom(p_local, n_local, prob_local,
		 lowertail_local, logp_local);
}

/* random generator */
void F77_SUB(rbinomrc)(double *n, double *prob, double *ans){
  double n_local = n[0];
  double prob_local = prob[0];
  ans[0] = rbinom(n_local, prob_local);
}

/* ------------------------------------------------------------------ */
/* poisson distribution */

/* mass function */
void F77_SUB(dpoisrc)(double *x, double *lambda, int *logd, double *ans){
  double x_local = x[0];
  double lambda_local = lambda[0];
  int logd_local = logd[0];
  ans[0] = dpois(x_local, lambda_local, logd_local);
}

/* mass function continuous version */
void F77_SUB(dpoisrawrc)(double *x, double *lambda, int *logd, double *ans){
  double x_local = x[0];
  double lambda_local = lambda[0];
  int logd_local = logd[0];
  ans[0] = dpois_raw(x_local, lambda_local, logd_local);
}

/* distribution function */
void F77_SUB(ppoisrc)(double *q, double *lambda,
		      int *lowertail, int *logp,  double *ans){
  double q_local = q[0];
  double lambda_local = lambda[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = ppois(q_local, lambda_local, lowertail_local, logp_local);
}

/* quantile function */
void F77_SUB(qpoisrc)(double *p, double *lambda,
		      int *lowertail, int *logp,  double *ans){
  double p_local = p[0];
  double lambda_local = lambda[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = qpois(p_local, lambda_local, lowertail_local, logp_local);
}

/* random generator */
void F77_SUB(rpoisrc)(double *lambda, double *ans){
  double lambda_local = lambda[0];
  ans[0] = rpois(lambda_local);
}

/* ------------------------------------------------------------------ */
/* negative binomial distribution */

/* mass function */
void F77_SUB(dnbinomrc)(double *x, double *n, double *prob, 
		      int *logd, double *ans){
  double x_local = x[0];
  double n_local = n[0];
  double prob_local = prob[0];
  int logd_local = logd[0];
  ans[0] = dnbinom(x_local, n_local, prob_local, logd_local);
}

/* distribution function */
void F77_SUB(pnbinomrc)(double *q, double *n, double *prob, 
		      int *lowertail, int *logp,  double *ans){
  double q_local = q[0];
  double n_local = n[0];
  double prob_local = prob[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = pnbinom(q_local, n_local, prob_local,
		 lowertail_local, logp_local);
}

/* quantile function */
void F77_SUB(qnbinomrc)(double *p, double *n, double *prob, 
		      int *lowertail, int *logp,  double *ans){
  double p_local = p[0];
  double n_local = n[0];
  double prob_local = prob[0];
  int lowertail_local = lowertail[0];
  int logp_local = logp[0];
  ans[0] = qnbinom(p_local, n_local, prob_local,
		 lowertail_local, logp_local);
}

/* random generator */
void F77_SUB(rnbinomrc)(double *n, double *prob, double *ans){
  double n_local = n[0];
  double prob_local = prob[0];
  ans[0] = rnbinom(n_local, prob_local);
}

/* ------------------------------------------------------------------ */
/* gamma function */

/* log gamma */
void F77_SUB(lgammarc)(double *x, double *ans){
  double x_local = x[0];
  ans[0] = lgammafn(x_local);
}

/* first deriv of log gamma */
void F77_SUB(digammarc)(double *x, double *ans){
  double x_local = x[0];
  ans[0] = digamma(x_local);
}

/* second deriv of log gamma */
void F77_SUB(trigammarc)(double *x, double *ans){
  double x_local = x[0];
  ans[0] = trigamma(x_local);
}

/* ------------------------------------------------------------------ */
/* misc math functions */

/* log(1 + x) */
void F77_SUB(log1prc)(double *x, double *ans){
  double x_local = x[0];
  ans[0] = log1p(x_local);
}

/* log(1 + x) - x */
void F77_SUB(log1pmxrc)(double *x, double *ans){
  double x_local = x[0];
  ans[0] = log1pmx(x_local);
}

/* log(1 + exp(x)) */
void F77_SUB(log1pexprc)(double *x, double *ans){
  double x_local = x[0];
  ans[0] = log1pexp(x_local);
}

/* exp(x) - 1 */
void F77_SUB(expm1rc)(double *x, double *ans){
  double x_local = x[0];
  ans[0] = expm1(x_local);
}

/* log(gamma(x + 1) */
void F77_SUB(lgamma1prc)(double *x, double *ans){
  double x_local = x[0];
  ans[0] = lgamma1p(x_local);
}
/* ------------------------------------------------------------------ */
