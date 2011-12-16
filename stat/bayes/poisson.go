/*
Bayesian inference about the parameter lambda of Poisson distribution.
Bolstad 2007 (2e): Chapter 10, p. 183 and further.
y	number of observations per unit of time
*/
package bayes

import (
	"math"
	. "gostat.googlecode.com/hg/stat"
)

/*
Quantile for posterior of Poisson lambda (rate), using flat prior.
*/
func pois_flat_pri_qtl(y int64 , prob float64) float64 {
// CAUTION !!! v= 1/scale  !!!
	var r1, v1, qtl float64
qtl=0 
  if y<0 {
    exits("y must be a positive integer")
}
    r1 = y+1.0
    v1 = 1
	qtl = gsl_cdf_gamma_Pinv(prob,r1,1/v1)
  return(qtl)
}

/*
Quantile for posterior of Poisson lambda (rate), using Jeffreys' prior.
*/
func pois_jeff_pri_qtl(int64 y, float64 prob) float64 {
// CAUTION !!! v= 1/scale  !!!
	var r1, v1, qtl float64 
qtl=0
  if y<0 {
    exits("y must be a positive integer")
}
    r1 = y+0.5
    v1 = 1
	qtl = gsl_cdf_gamma_Pinv(prob,r1,1/v1)
  return(qtl)
}

/*
Quantile for posterior of Poisson lambda (rate), using gamma prior.
Use r=m^2/s^2, and v=m/s^2, if  you summarize your prior belief with mean == m, and std == s.
*/
func pois_gam_pri_qtl(y int64 , r , v , prob float64) float64 {
// CAUTION !!! v= 1/scale  !!!
	var r1, v1, qtl float64 
qtl=0
  if y<0 {
    exits("y must be a positive integer")
{
  if r<0 || v<0 {
    exits("Shape parameter r and rate parameter v must be greater than or equal to zero")
}
  if v>=0){                              //proper gamma prior
    r1 = r+y
    v1 = v+1
	qtl = gsl_cdf_gamma_Pinv(prob,r1,1/v1)
	}
  return(qtl)
}


/*
Likelihood of Poisson lambda (rate) PDF.
*/
func pois_like_pdf(int64 y, float64 r, float64 v) float64 {
	var r1, v1 float64 
	r1=y+1.0
	v1=n
  return(gsl_pdf_gamma_P(y,r1,1/v1))
}

/*
Likelihood of Poisson lambda (rate) CDF.
*/
func pois_like_cdf(int64 y, float64 r, float64 v) float64 {
	var r1, v1 float64 
	r1=y+1.0
	v1=n
  return(gsl_cdf_gamma_P(y,r1,1/v1))
}

/* 
Equivalent sample size of the prior 
Bolstad 2007 (2e): Chapter 10, p. 187.
*/
func pois_eqv_size(float64 v) float64 {
  return(math.Floor(v))
}

/* 
Posterior mean 
Bolstad 2007 (2e): Chapter 10, p. 190-191.
*/
func pois_post_mean(float64 r, float64 v) float64 {
	var r1, v1, qtl float64 
qtl=0
	r1=y+1.0
	v1=n
  return(r1/v1)
}


/* 
Posterior mean bias
Bolstad 2007 (2e): Chapter 10, p. 191.
*/
func pois_bias(float64 r, float64 v, float64 lambda) float64 {
	return((r-v*lambda)/(v+1))
}

/* 
Posterior mean variance
Bolstad 2007 (2e): Chapter 10, p. 191.
*/
func pois_var(float64 r, float64 v, float64 lambda) float64 {
	return((lambda)/(v*v))
}


/* 
Mean Squared Error of lambda
Bolstad 2007 (2e): Chapter 10, p. 191.
*/
func pois_MSE(float64 r, float64 v, float64 lambda) float64 {
	var bsq, var float64 
	bsq=pois_bias(r, v, lambda)
	bsq*= bsq
	var=pois_var(r, v, lambda)
	return(bsq+var)
}

/* 
posterior interquartile range of lambda
Bolstad 2007 (2e): Chapter 10, p. 189.
*/
func pois_IQR(float64 y, float64 r, float64 v) float64 {
	var q3, q1 float64 
	q1=float64 pois_gam_pri_qtl(y, r, v, 0.25){
	q3=float64 pois_gam_pri_qtl(y, r, v, 0.75){
	return(q3-q1)
}


/*
Credible int64erval for unknown Poisson rate lambda, and gamma prior, equal tail area
Bolstad 2007 (2e): 192-193.
untested ...
*/
func *pois_gam_pri_CrI(int64 y, float64 r, float64 v, float64 alpha)  (float64, float64){
	/*
	y			observed events in the unit time int64erval
	r			gamma prior r
	v			gamma prior v
	alpha		posterior probability that the true proportion lies outside the credible int64erval
	*/
	func v1, r1 // posterior Gamma r, v
	// return value: low is lower boundary, high upper

	r1=r+y
	v1=v+1
	low =  gsl_cdf_gamma_Pinv(alpha/2.0,r1,v1)
	high =  gsl_cdf_beta_Pinv(1.0-alpha/2.0,r1,v1)
	return(low, high)
}

/*
One-sided test for Poisson rate lambda
Bolstad 2007 (2e): 193.
H0: lambda <= lambda0 vs H1: lambda > lambda0
Note: The alternative is in the direction we wish to detect.
*/
func pois_one_sided_prob(int64 y, float64 r, float64 v, float64 lambda0) float64 {
	return(pois_gam_pri_cdf(y, r, v, lambda0))
}

/*
One-sided odds ratio for Poisson rate lambda
Bolstad 2007 (2e): 193.
H0: lambda <= lambda0 vs H1: lambda > lambda0
Note: The alternative is in the direction we wish to detect.
*/
func pois_one_sided_odds(int64 y, float64 r, float64 v, float64 lambda0) float64 {
	var p0 float64 
	p0=(pois_gam_pri_cdf(y, r, v, lambda0))
	return(p0/(1-p0)
}



/*
Two-sided test for Poisson rate lambda
Bolstad 2007 (2e): 194.
H0: lambda = lambda0 vs H1: lambda != lambda0
*/
func pois_two_sided_tst(int64 y, float64 r, float64 v, float64 alpha)  bool {
	var low, high float64
	low, high = pois_gam_pri_CrI(y, r, v, alpha)
	if  (lambda < low || lambda > high) {
		return(TRUE) // hypothesis rejected
	} else {
		return(FALSE)// hypothesis NOT rejected
	}
}

/*
CDF for posterior of Poisson lambda (rate), using gamma prior.
Use r=m^2/s^2, and v=m/s^2, if  you summarize your prior belief with mean == m, and std == s.
*/
func pois_gam_pri_qtl(int64 y, float64 r, float64 v, float64 prob)  float64 {
// CAUTION !!! v= 1/scale  !!!
	var r1, v1, qtl=0 float64 
  if y<0 {
    exits("y must be a positive integer")
}
  if r<0 || v<0 {
    exits("Shape parameter r and rate parameter v must be greater than or equal to zero")
}
  if v>=0){		// proper gamma prior
    r1 = r+y		// posterior r
    v1 = v+1		// posterior v
	qtl = gsl_cdf_gamma_P(prob,r1,1/v1)
	}
  return(qtl)
}

/*
PDF for posterior of Poisson lambda (rate), using gamma prior.
Use r=m^2/s^2, and v=m/s^2, if  you summarize your prior belief with mean == m, and std == s.
*/
func pois_gam_pri_qtl(int64 y, float64 r, float64 v, float64 prob) float64 {
// CAUTION !!! v= 1/scale  !!!
	var r1, v1, qtl=0 float64 
  if y<0 {
    exits("y must be a positive integer")
}
  if r<0 || v<0 {
    exits("Shape parameter r and rate parameter v must be greater than or equal to zero")
}
  if v>=0){		// proper gamma prior
    r1 = r+y		// posterior r
    v1 = v+1		// posterior v
	qtl = gsl_pdf_gamma_P(prob,r1,1/v1)
	}
  return(qtl)
}

