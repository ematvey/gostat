// Bayesian inference about the parameter λ of Poisson distribution.
// Bolstad 2007 (2e): Chapter 10, p. 183 and further.
// y	number of observations per unit of time

package bayes

import (
	. "gostat.googlecode.com/hg/stat"
	"math"
)

//Quantile for posterior of Poisson λ (rate), using flat prior.
func PoissonLambda_Qtl_FPri(y int64, prob float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	var r1, v1 float64
	if y < 0 {
		panic("y must be a positive integer")
	}
	r1 = float64(y) + 1.0
	v1 = 1
	return Gamma_Qtl_For(r1, 1/v1, prob)
}

// Quantile for posterior of Poisson λ (rate), using Jeffreys' prior.
func PoissonLambda_Qtl_JPri(y int64, prob float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	var r1, v1 float64
	if y < 0 {
		panic("y must be a positive integer")
	}
	r1 = float64(y) + 0.5
	v1 = 1
	return Gamma_Qtl_For(r1, 1/v1, prob)
}

// Quantile for posterior of Poisson λ (rate), using gamma prior.
// Use r=m^2/s^2, and v=m/s^2, if you summarize your prior belief with mean == m, and std == s.
func PoissonLambda_Qtl_GPri(y int64, r, v, prob float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	var r1, v1 float64
	if y < 0 {
		panic("y must be a positive integer")
	}
	if r < 0 || v < 0 {
		panic("Shape parameter r and rate parameter v must be greater than or equal to zero")
	}
	r1 = r + float64(y)
	v1 = v + 1
	return Gamma_Qtl_For(r1, 1/v1, prob)
}

/*
// Likelihood of Poisson λ (rate) PDF.
func PoissonLambda_Like_PDF(y, n int64, r, v float64) float64 {
	var r1, v1 float64
	r1 = float64(y) + 1.0
	v1 = float64(n)
	return Gamma_PDF_At(r1, 1/v1, float64(y))
	//poisson.go:60: internal compiler error: unknown type
}

// Likelihood of Poisson λ (rate) CDF.
func PoissonLambda_Like_CDF(y, n int64, r, v float64) float64 {
	var r1, v1 float64
	r1 = float64(y) + 1.0
	v1 = float64(n)
	return Gamma_CDF_At(r1, 1/v1, float64(y))
}
*/


// Equivalent sample size of the prior 
// Bolstad 2007 (2e): Chapter 10, p. 187.
func PoissonLambdaEqvSize(v float64) float64 {
	return (math.Floor(v))
}

// Posterior mean 
// Bolstad 2007 (2e): Chapter 10, p. 190-191.
func PoissonLambdaPostMean(y, n int64, r, v float64) float64 {
	var r1, v1 float64
	r1 = float64(y) + 1.0
	v1 = float64(n)
	return r1 / v1
}

// Posterior mean bias
// Bolstad 2007 (2e): Chapter 10, p. 191.
func PoissonLambdaPostMeanBias(r, v, λ float64) float64 {
	return (r - v*λ) / (v + 1)
}

// Posterior variance
// Bolstad 2007 (2e): Chapter 10, p. 191.
func PoissonLambdaPostVar(r, v, λ float64) float64 {
	return λ / (v * v)
}

// Mean Squared Error of λ
// Bolstad 2007 (2e): Chapter 10, p. 191.
func PoissonLambdaMSE(r, v, λ float64) float64 {
	var bsq, variance float64
	bsq = PoissonLambdaPostMeanBias(r, v, λ)
	bsq *= bsq
	variance = PoissonLambdaPostVar(r, v, λ)
	return (bsq + variance)
}

// posterior interquartile range of λ
// Bolstad 2007 (2e): Chapter 10, p. 189.
func PoissonLambdaIQR(y int64, r, v float64) float64 {
	var q3, q1 float64
	q1 = PoissonLambda_Qtl_GPri(y, r, v, 0.25)
	q3 = PoissonLambda_Qtl_GPri(y, r, v, 0.75)
	return (q3 - q1)
}

// Credible interval for unknown Poisson rate λ, and gamma prior, equal tail area
// Bolstad 2007 (2e): 192-193.
// untested ...
func PoissonLambda_CrI_GPri(y int64, r, v, α float64) (float64, float64) {
	/*
		y			observed events in the unit time interval
		r			gamma prior r
		v			gamma prior v
		α		posterior probability that the true proportion lies outside the credible interval
	*/
	// return value: low is lower boundary, high upper

	r1 := r + float64(y) // posterior Gamma r
	v1 := v + 1          // posterior Gamma v
	low := Gamma_Qtl_For(r1, v1, α/2.0)
	high := Gamma_Qtl_For(r1, v1, 1.0-α/2.0)
	return low, high
}

/*
One-sided test for Poisson rate λ
Bolstad 2007 (2e): 193.
H0: λ <= λ0 vs H1: λ > λ0
Note: The alternative is in the direction we wish to detect.
*/
func PoissonLambda_OneSidedP(y int64, r, v, λ0 float64) float64 {
	return pois_gam_pri_cdf(y, r, v, λ0)
}

/*
One-sided odds ratio for Poisson rate λ
Bolstad 2007 (2e): 193.
H0: λ <= λ0 vs H1: λ > λ0
Note: The alternative is in the direction we wish to detect.
*/
func PoissonLambda_OneSidedOdds(y int64, r, v, λ0 float64) float64 {
	var p0 float64
	p0 = pois_gam_pri_cdf(y, r, v, λ0)
	return p0 / (1 - p0)
}

/*
Two-sided test for Poisson rate λ
Bolstad 2007 (2e): 194.
H0: λ = λ0 vs H1: λ != λ0
func PoissonLambda_TwoSidedTst(y int64, r, v, α, λ float64) bool {
	var low, high float64
	low, high = PoissonLambda_CrI_GPri(y, r, v, α)
	if λ < low || λ > high {
		return (TRUE) // hypothesis rejected
	} else {
		return (FALSE) // hypothesis NOT rejected
	}
}
*/

/*
CDF for posterior of Poisson λ (rate), using gamma prior.
Use r=m^2/s^2, and v=m/s^2, if you summarize your prior belief with mean == m, and std == s.
*/
func pois_gam_pri_cdf(y int64, r, v, prob float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	var r1, v1 float64
	if y < 0 {
		panic("y must be a positive integer")
	}
	if r < 0 || v < 0 {
		panic("Shape parameter r and rate parameter v must be greater than or equal to zero")
	}
	r1 = r + float64(y) // posterior r
	v1 = v + 1          // posterior v
	return Gamma_CDF_At(r1, 1/v1, prob)
}

/*
PDF for posterior of Poisson λ (rate), using gamma prior.
Use r=m^2/s^2, and v=m/s^2, if you summarize your prior belief with mean == m, and std == s.
*/
func pois_gam_pri_pdf(y int64, r, v, prob float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	var r1, v1 float64
	if y < 0 {
		panic("y must be a positive integer")
	}
	if r < 0 || v < 0 {
		panic("Shape parameter r and rate parameter v must be greater than or equal to zero")
	}
	r1 = r + float64(y) // posterior r
	v1 = v + 1          // posterior v
	return Gamma_PDF_At(r1, 1/v1, prob)
}
