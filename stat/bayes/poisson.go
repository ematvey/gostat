// Bayesian inference about the parameter λ of Poisson distribution.
// Bolstad 2007 (2e): Chapter 10, p. 183 and further.
// sumK	sum of observations over all repetitions
// n number of repetitions

package bayes

import (
	. "gostat.googlecode.com/hg/stat/prob"
//	. "go-fn.googlecode.com/hg/fn"
	"math"
)

// Poisson λ, posterior PDF, flat prior.
func PoissonLambda_PDF_FPri(sumK, n int64) func(p float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	if sumK < 0 || n <= 0 {
		panic("bad data")
	}
	r1 := float64(sumK) + 1.0
	v1 := float64(n)
	return Gamma_PDF(r1, 1/v1)
}
// Poisson λ, posterior PDF, Jeffreys' prior.
func PoissonLambda_PDF_JPri(sumK, n int64) func(p float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	if sumK < 0 || n <= 0 {
		panic("bad data")
	}
	r1 := float64(sumK) + 0.5
	v1 := float64(n)
	return Gamma_PDF(r1, 1/v1)
}

// Poisson λ, posterior PDF, gamma prior.
// Use r=m^2/s^2, and v=m/s^2, if you summarize your prior belief with mean == m, and std == s.
func PoissonLambda_PDF_GPri(sumK, n int64, r, v float64) func(p float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	if sumK < 0 || n <= 0 {
		panic("bad data")
	}
	if r < 0 || v < 0 {
		panic("Shape parameter r and rate parameter v must be greater than or equal to zero")
	}
	r1 := r + float64(sumK)
	v1 := v + float64(n)
	return Gamma_PDF(r1, 1/v1)
}

// Poisson λ, posterior CDF, flat prior.
func PoissonLambda_CDF_FPri(sumK, n int64) func(p float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	if sumK < 0 || n <= 0 {
		panic("bad data")
	}
	r1 := float64(sumK) + 1.0
	v1 := float64(n)
	return Gamma_CDF(r1, 1/v1)
}
// Poisson λ, posterior CDF, Jeffreys' prior.
func PoissonLambda_CDF_JPri(sumK, n int64) func(p float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	if sumK < 0 || n <= 0 {
		panic("bad data")
	}
	r1 := float64(sumK) + 0.5
	v1 := float64(n)
	return Gamma_CDF(r1, 1/v1)
}

// Poisson λ, posterior CDF, gamma prior.
// Use r=m^2/s^2, and v=m/s^2, if you summarize your prior belief with mean == m, and std == s.
func PoissonLambda_CDF_GPri(sumK, n int64, r, v float64) func(p float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	if sumK < 0 || n <= 0 {
		panic("bad data")
	}
	if r < 0 || v < 0 {
		panic("Shape parameter r and rate parameter v must be greater than or equal to zero")
	}
	r1 := r + float64(sumK)
	v1 := v + float64(n)
	return Gamma_CDF(r1, 1/v1)
}

// Poisson λ, posterior quantile function, flat prior.
func PoissonLambda_Qtl_FPri(sumK, n int64) func(p float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	if sumK < 0 || n <= 0 {
		panic("bad data")
	}
	r1 := float64(sumK) + 1.0
	v1 := float64(n)
	return Gamma_Qtl(r1, 1/v1)
}
// Poisson λ, posterior quantile function, Jeffreys' prior.
func PoissonLambda_Qtl_JPri(sumK, n int64) func(p float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	if sumK < 0 || n <= 0 {
		panic("bad data")
	}
	r1 := float64(sumK) + 0.5
	v1 := float64(n)
	return Gamma_Qtl(r1, 1/v1)
}

// Poisson λ, posterior quantile function, gamma prior.
// Use r=m^2/s^2, and v=m/s^2, if you summarize your prior belief with mean == m, and std == s.
func PoissonLambda_Qtl_GPri(sumK, n int64, r, v float64) func(p float64) float64 {
	// CAUTION !!! v= 1/scale !!!
	if sumK < 0 || n <= 0 {
		panic("bad data")
	}
	if r < 0 || v < 0 {
		panic("Shape parameter r and rate parameter v must be greater than or equal to zero")
	}
	r1 := r + float64(sumK)
	v1 := v + float64(n)
	return Gamma_Qtl(r1, 1/v1)
}

// Likelihood of Poisson λ.
// Bolstad 2007 (2e): Chapter 10, p. 184.
func PoissonLambda_Like(sumK, n int64, λ float64) float64 {
	return λ*float64(sumK)* math.Exp(float64(-n)*λ)

}

// Equivalent sample size of the prior 
// Bolstad 2007 (2e): Chapter 10, p. 187.
func PoissonLambdaEqvSize(v float64) float64 {
	return (math.Floor(v))
}

// Posterior mean 
// Bolstad 2007 (2e): Chapter 10, p. 190-191.
func PoissonLambdaPostMean(sumK, n int64, r, v float64) float64 {
	r1 := float64(sumK) + 1.0
	v1 := 1.0
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
	bsq := PoissonLambdaPostMeanBias(r, v, λ)
	bsq *= bsq
	variance := PoissonLambdaPostVar(r, v, λ)
	return (bsq + variance)
}

// posterior interquartile range of λ
// Bolstad 2007 (2e): Chapter 10, p. 189.
func PoissonLambdaIQR(sumK, n int64, r, v float64) float64 {
	qf := PoissonLambda_Qtl_GPri(sumK, n, r, v)
	q1 := qf(0.25)
	q3 := qf(0.75)
	return (q3 - q1)
}

// Credible interval for unknown Poisson rate λ, and gamma prior, equal tail area
// Bolstad 2007 (2e): 192-193.
// untested ...
func PoissonLambda_CrI_GPri(sumK, n int64, r, v, α float64) (lo, hi float64) {
	/*
		sumK, n			total observed events in n equal time intervals
		r			gamma prior r
		v			gamma prior v
		α		posterior probability that the true proportion lies outside the credible interval
	*/
	// return value: lo is lower boundary, hi upper
	qf := PoissonLambda_Qtl_GPri(sumK, n, r, v)
	lo = qf(α/2)
	hi = qf(1-α/2)
	return
}
// One-sided test for Poisson rate λ
// Bolstad 2007 (2e): 193.
// H0: λ <= λ0 vs H1: λ > λ0
// Note: The alternative is in the direction we wish to detect.
func PoissonLambda_OneSidedTst(sumK, n int64, r, v, α, λ0 float64) bool {
	cdf := PoissonLambda_CDF_GPri(sumK, n, r, v)
	p0 := cdf(λ0)
	reject := false // hypothesis NOT rejected (default)
	if p0 < α {
		reject = true // hypothesis rejected
	}
	return reject
}

// One-sided odds ratio for Poisson rate λ
// Bolstad 2007 (2e): 193.
// H0: λ <= λ0 vs H1: λ > λ0
// Note: The alternative is in the direction we wish to detect.
func PoissonLambda_OneSidedOdds(sumK, n int64, r, v, λ0 float64) float64 {
	cdf := PoissonLambda_CDF_GPri(sumK, n, r, v)
	p0 := cdf(λ0)
	return p0 / (1 - p0)
}

// Two-sided test for Poisson rate λ
// Bolstad 2007 (2e): 194.
// H0: λ = λ0 vs H1: λ != λ0
func PoissonLambda_TwoSidedTst(sumK, n int64, r, v, α, λ0 float64) bool {
	low, high := PoissonLambda_CrI_GPri(sumK, n, r, v, α)
	reject := false // hypothesis NOT rejected (default)
	if λ0 < low || λ0 > high {
		reject = true // hypothesis rejected
	}
	return reject
}


