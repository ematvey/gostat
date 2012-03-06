// Bayesian inference about the parameter p of binomial distribution.
// Bolstad 2007 (2e): Chapter 8, p. 141 and further.

package bayes

import (
	"fmt"
	. "gostat.googlecode.com/hg/stat/prob"
	"math"
)

// Quantile, Flat prior 
func BinomPi_Qtl_FPri(k, n int64, p float64) float64 {
	var α, β float64
	if k > n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}

	α = 1.0
	β = 1.0
	return Beta_Qtl_For(α+float64(k), β+float64(n-k), p)
}

// Quantile, Jeffreys prior
// see Aitkin 2010: 143 for cautions
func BinomPi_Qtl_JPri(k, n int64, p float64) float64 {
	var α, β float64
	α = 0.5
	β = 0.5
	if k > n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}
	return Beta_Qtl_For(α+float64(k), β+float64(n-k), p)
}

// Quantile, Haldane prior
// see Aitkin 2010: 143 for cautions
func BinomPi_Qtl_HPri(k, n int64, p float64) float64 {
	var α, β float64
	α = 0.0
	β = 0.0
	if k > n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}
	return Beta_Qtl_For(α+float64(k), β+float64(n-k), p)
}

// Quantile, Beta prior, general
func BinomPi_Qtl_BPri(k, n int64, α, β, p float64) float64 {
	if k > n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}
	if α < 0 || β < 0 {
		panic(fmt.Sprintf("The parameters of the prior must be non-negative"))
	}
	return Beta_Qtl_For(α+float64(k), β+float64(n-k), p)
}

// Equivalent sample size of the prior
func BinomPiEqvSize(α, β float64) int64 {
	return int64(math.Floor(α + β + 1))
}

// Posterior modus
func BinomPiPostModus(α, β float64, n, k int64) float64 {
	var post_α, post_β float64
	post_α = α + float64(k)
	post_β = β + float64(n-k)
	return (post_α - 1) / (post_α + post_β - 2.0)
}

//  Posterior mean
func BinomPiPostMean(α, β float64, n, k int64) float64 {
	var post_α, post_β float64
	post_α = α + float64(k)
	post_β = β + float64(n-k)
	return ((post_α) / (post_α + post_β))
}

//  Posterior median
func BinomPiPostMedian(α, β float64, n, k int64) float64 {
	return 0 // to be implemented
}

// Posterior variance
// Bolstad 2007 (2e): 151, eq. 8.5
func BinomPiPostVar(α, β float64, n, k int64) float64 {
	var post_α, post_β float64
	post_α = α + float64(k)
	post_β = β + float64(n-k)
	return (post_α * post_β) / ((post_α + post_β) * (post_α + post_β) * (post_α + post_β + 1.0))
}

//  Posterior mean square of p
// Bolstad 2007 (2e): 152-153, eq. 8.7
func BinomPiPMS(α, β float64, n, k, which_pi int64) float64 {
	const (
		MEAN   = 0
		MEDIAN = 1
		MODUS  = 2
	)

	var post_mean, post_var, pi_hat float64

	post_var = BinomPiPostVar(α, β, n, k)
	post_mean = BinomPiPostMean(α, β, n, k)

	switch which_pi {
	case MEAN:
		pi_hat = BinomPiPostMean(α, β, n, k)
	case MEDIAN:
		pi_hat = BinomPiPostMedian(α, β, n, k)
	case MODUS:
		pi_hat = BinomPiPostModus(α, β, n, k)
	}
	return post_var + (post_mean-pi_hat)*(post_mean-pi_hat)
}

// Credible interval for unknown binomial proportion, and beta prior, equal tail area
// Bolstad 2007 (2e): 153
// untested ...

func BinomPi_CrI_BP(α, β, alpha float64, n, k int64) (float64, float64) {
	/*
		k			observed successes
		n			total number of observations
		α			beta prior a
		β			beta prior b
		alpha			posterior probability that the true proportion lies outside the credible interval
	*/

	var low, upp float64
	low = Beta_Qtl_For(alpha/2.0, α+float64(k), β+float64(n-k))
	upp = Beta_Qtl_For(1.0-alpha/2.0, α+float64(k), β+float64(n-k))
	return low, upp
}

// Credible interval for unknown binomial proportion, and beta prior, equal tail area, normal approximation
// Bolstad 2007 (2e): 154-155, eq. 8.8
// untested ...

func BinomPi_CrI_BPriNApprox(α, β, alpha float64, n, k int64) (float64, float64) {
	/*
		k			observed successes
		n			total number of observations
		a			beta prior a
		b			beta prior b
		alpha			posterior probability that the true proportion lies outside the credible interval
	*/

	var post_mean, post_var, post_α, post_β, z, low, upp float64

	post_α = α + float64(k)
	post_β = β + float64(n-k)

	post_mean = post_α / (post_α + post_β)
	post_var = (post_α * post_β) / ((post_α + post_β) * (post_α + post_β) * (post_α + post_β + 1.0))
	z = Z_Qtl_For(alpha / 2)

	low = post_mean - z*math.Sqrt(post_var)
	upp = post_mean + z*math.Sqrt(post_var)
	return low, upp
}


// Likelihood
func BinomPi_Like(pi float64, n, k int64) float64 {
	return math.Pow(pi, float64(k)) * math.Pow(1-pi, float64(n-k))
}

// Deviance 
func BinomPi_Deviance(pi float64, n, k int64) float64 {
	return -2*math.Log(BinomPi_Like(pi, n, k))
}

// Posterior PDF, Beta prior
func BinomPi_PDF_BPri(k, n int64, α, β, p float64) float64 {
	if k > n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}
	if α < 0 || β < 0 {
		panic(fmt.Sprintf("The parameters of the prior must be non-negative"))
	}
	return Beta_PDF_At(α+float64(k), β+float64(n-k), p)
}

// Posterior CDF, Beta prior
func BinomPi_CDF_BPri(k, n int64, α, β, p float64) float64 {
	if k > n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}
	if α < 0 || β < 0 {
		panic(fmt.Sprintf("The parameters of the prior must be non-negative"))
	}
	return Beta_CDF_At(α+float64(k), β+float64(n-k), p)
}

// Sampling from posterior, Beta prior
func NextBinomPi_CDF_BPri(k, n int64, α, β float64) float64 {
	if k > n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}
	if α < 0 || β < 0 {
		panic(fmt.Sprintf("The parameters of the prior must be non-negative"))
	}
	return NextBeta(α+float64(k), β+float64(n-k))
}

// Deviance difference of a point null hypothesis pi = p against general alternative pi != p
// Aitkin 2010:143-144.
func binomPi_PointDevDiff(k, n int64, α, β, p, pi float64) float64 {
	nn := float64(n)
	kk := float64(k)
	d0 := -2*(kk*math.Log(p)+(nn-kk)*math.Log(1-p))	//  null model deviance
	dd := d0 + 2*(kk*math.Log(pi)+(nn-kk)*math.Log(1-pi))
	return dd
}


	
