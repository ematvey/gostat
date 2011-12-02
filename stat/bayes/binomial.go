// Bayesian inference about the parameter p of binomial distribution.
// Bolstad 2007 (2e): Chapter 8, p. 141 and further.

package bayes
import (
	"fmt"
	"math"
	s "gostat.googlecode.com/hg/stat"
)


// Quantile, Flat prior 
func BinomFlatPriQtl(k, n int64, p float64) float64 {
	var α, β float64
	if k>n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}

	α=float64(k+1)
	β=float64(n-k+1)

	if α<=0||β<=0 {
		panic(fmt.Sprintf("The parameters of the prior must be greater than zero"))
	}
	return	s.BetaInv_CDF_For(α, β, p)
}

// Quantile, Beta prior
func BinomBetaPriQtl(k, n int64, α, β, p float64) float64 {
	if k>n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}
	if α<=0||β<=0 {
		panic(fmt.Sprintf("The parameters of the prior must be greater than zero"))
	}
	return	s.BetaInv_CDF_For(α+float64(k), β+float64(n-k), p)
}

// Quantile, Jeffrey's prior
func BinomJeffPriQtl(k, n int64, p float64) float64 {
	var α, β float64
	α=0.5
	β=0.5
	if k>n {
		panic(fmt.Sprintf("The number of observed successes (k) must be <= number of trials (n)"))
	}
	return	s.BetaInv_CDF_For(α+float64(k), β+float64(n-k), p)
}

// Equivalent sample size of the prior
func BinomEqvSize(α,  β float64) int64 {
	return int64(math.Floor(α+β+1))
}

// Posterior modus
func BinomPostModus(α,  β float64, n, k int64) float64 {
	var post_α, post_β float64
	post_α=α+float64(k)
	post_β=β+float64(n-k)
	return (post_α-1)/(post_α+post_β-2.0)
}


//  Posterior mean
func BinomPostMean(α,  β float64, n, k int64) float64 {
	var post_α, post_β float64
	post_α=α+float64(k)
	post_β=β+float64(n-k)
	return((post_α)/(post_α+post_β))
}

//  Posterior median
func BinomPostMedian(α,  β float64, n, k int64) float64 {
	return 0  // to be implemented
}

// Posterior variance
// Bolstad 2007 (2e): 151, eq. 8.5
func BinomPostVar(α, β float64, n, k int64) float64 {
	var post_α, post_β float64
	post_α=α+float64(k)
	post_β=β+float64(n-k)
	return (post_α*post_β)/((post_α+post_β)*(post_α+post_β)*(post_α+post_β+1.0))
}

//  Posterior mean square of p
// Bolstad 2007 (2e): 152-153, eq. 8.7
func BinomPMS(α,  β float64, n, k, which_pi int64) float64 {
	const (
		MEAN = 0
		MEDIAN = 1
		MODUS = 2
	)

	var post_mean, post_var, pi_hat float64

	post_var=BinomPostVar(α, β , n, k )
	post_mean = BinomPostMean(α, β , n, k )

	switch(which_pi) {
	case MEAN : 		pi_hat=BinomPostMean(α,  β , n, k)
	case MEDIAN : 		pi_hat=BinomPostMedian(α,  β , n, k)
	case MODUS : 		pi_hat=BinomPostModus(α,  β , n, k)
	}
	return post_var+ (post_mean-pi_hat)*(post_mean-pi_hat)
}

// Credible interval for unknown binomial proportion, and beta prior, equal tail area
// Bolstad 2007 (2e): 153
// untested ...


func BinomBetaPriCrI(α,  β, alpha float64, n, k int64)  (float64, float64) {
	/*
	k			observed successes
	n			total number of observations
	α			beta prior a
	β			beta prior b
	alpha			posterior probability that the true proportion lies outside the credible interval
	*/

	var low, upp float64
	low =  s.BetaInv_CDF_For(alpha/2.0,α+float64(k),β+float64(n-k))
	upp =  s.BetaInv_CDF_For(1.0-alpha/2.0,α+float64(k),β+float64(n-k))
	return low, upp
}

// Credible interval for unknown binomial proportion, and beta prior, equal tail area, normal approximation
// Bolstad 2007 (2e): 154-155, eq. 8.8
// untested ...

func BinomBetaPriCrIApprox(α,  β, alpha float64, n, k int64) (float64, float64) {
	/*
	k			observed successes
	n			total number of observations
	a			beta prior a
	b			beta prior b
	alpha			posterior probability that the true proportion lies outside the credible interval
	*/

	var post_mean, post_var, post_α, post_β, z, low, upp float64


	post_α=α+float64(k)
	post_β=β+float64(n-k)

	post_mean = post_α/(post_α+post_β)
	post_var = (post_α*post_β)/((post_α+post_β)*(post_α+post_β)*(post_α+post_β+1.0))
	z = s.ZInv_CDF_For(alpha/2)

	low =  post_mean - z * math.Sqrt(post_var)
	upp =  post_mean + z * math.Sqrt(post_var)
	return low, upp
}


