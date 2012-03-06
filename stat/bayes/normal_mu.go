// Bayesian inference about the params (mu, sigma) of Normal (Gaussian) distribution.
// Bolstad 2007 (2e): Chapter 11, p. 199 and further.

package bayes

import (
	"fmt"
	. "gostat.googlecode.com/hg/stat/prob"
	"math"
)

// PMF of the posterior distribution of unknown Normal mean, with KNOWN sigma, and discrete prior, for single observation
// Bolstad 2007 (2e): 200-201
// y		single observation taken from Normal distribution
// sigma	standard deviation of population, assumed to be known
// p		probability for which the quantile will be returned
func NormMuSingle_PMF_DPri(y, sigma float64, mu []float64, prior []float64) (post []float64) {
	n := len(mu)
	if len(prior) != n {
		panic(fmt.Sprintf("Discrete means and their priors must be vectors of the same length"))
	}
	post = make([]float64, n)
	sum := 0.0
	for i := 0; i< n; i++ {
		z := (y-mu[i])/sigma
		like := Z_PDF_At(z)
		post[i] = prior[i] * like
		sum += post[i]
	}
	for i := 0; i< n; i++ {
		post[i] /= sum
	}
	return
}

// Posterior mean for unknown Normal mean, with KNOWN sigma. 
// Bolstad 2007 (2e): 209, eq. 11.6
func NormMuPostMean(n int, samp_mean, sigma, pri_mu, pri_sigma float64) float64 {
	// pri_mu		prior mean
	// pri_sigma		prior standard deviation
	// n			size of sample == number of measurements
	// sigma		standard deviation of population, assumed to be known (alternatively, use an estimate)
	σ2 := sigma * sigma
	nn := float64(n)
	pri_var := pri_sigma * pri_sigma
	post_mu := (pri_mu/pri_var)/(nn/σ2+1/pri_var) + samp_mean*(nn/σ2)/(nn/σ2+1/pri_var)
	return (post_mu)
}

// Posterior standard deviation for unknown Normal mean, with KNOWN sigma. 
// Bolstad 2007 (2e): 209, eq. 11.6
func NormMuPostStd(n int, sigma, pri_mu, pri_sigma float64) float64 {
	// pri_mu		prior mean
	// pri_sigma		prior standard deviation
	// n			size of sample == number of measurements
	// sigma		standard deviation of population, assumed to be known (alternatively, use an estimate)
	σ2 := sigma * sigma
	nn := float64(n)
	pri_var := pri_sigma * pri_sigma
	post_var := (σ2 * pri_var) / (σ2 + float64(nn)*pri_var)
	post_sigma := math.Sqrt(post_var)
	return (post_sigma)
}

// Quantile for posterior distribution of unknown Normal mean, with KNOWN sigma, and flat prior (Jeffrey's prior), for single observation
// Bolstad 2007 (2e): 206
// y		single observation taken from Normal distribution
// sigma	standard deviation of population, assumed to be known
// p		probability for which the quantile will be returned
// untested ...
func NormMuSingle_Qtl_FPri(y, sigma, p float64) float64 {
	post_mu := y
	post_sigma := sigma
	qtl := Normal_Qtl_For(post_mu, post_sigma, p)
	return (qtl)
}

// Quantile for posterior distribution of unknown Normal mean, with KNOWN sigma, and flat prior (Jeffrey's prior), for sample
// Bolstad 2007 (2e): 207
// samp_mean		sample mean of observations taken from Normal distribution
// sigma		standard deviation of population, assumed to be known
// n			number of observations
// p			probability for which the quantile will be returned
func NormMu_Qtl_FPri(n int, samp_mean, sigma, p float64) float64 {
	// check params
	if sigma <= 0 {
		panic(fmt.Sprintf("Prior standard deviation must be greater than zero"))
	}

	nn := float64(n)
	σ2 := sigma * sigma
	post_mu := samp_mean
	post_var := σ2 / nn
	post_sigma := math.Sqrt(post_var)
	return Normal_Qtl_For(post_mu, post_sigma, p)
}

// Quantile for posterior distribution of unknown Normal mean, with KNOWN sigma, and normal prior, for single observation
// Bolstad 2007 (2e): 208, eq. 11.4
// y		single observation taken from Normal distribution
// sigma	standard deviation of population, assumed to be known
// pri_mu	Normal prior mean
// pri_sigma	Normal prior standard deviation
// p		probability for which the quantile will be returned
// untested ...
func NormMuSingle_Qtl_NPri(y, sigma, pri_mu, pri_sigma, p float64) float64 {
	σ2 := sigma * sigma
	pri_var := pri_sigma * pri_sigma
	post_mu := (σ2*pri_mu + pri_var*y) / (σ2 + pri_var)
	post_var := (σ2 * pri_var) / (σ2 + pri_var)
	post_sigma := math.Sqrt(post_var)
	return Normal_Qtl_For(post_mu, post_sigma, p)
}

// Quantile for posterior distribution of unknown Normal mean, with KNOWN sigma, and normal prior, for sample
// Bolstad 2007 (2e): 209, eq. 11.5, 11.6
// samp_mean		sample mean of observations taken from Normal dist
// sigma		standard deviation of population, assumed to be known
// n			number of observations
// pri_mu		Normal prior mean
// pri_sigma		Normal prior standard deviation
// p			probability for which the quantile will be returned
func NormMu_Qtl_NPri(n int, samp_mean, sigma, pri_mu, pri_sigma, p float64) float64 {
	nn := float64(n)
	σ2 := sigma * sigma
	pri_var := pri_sigma * pri_sigma
	post_var := (σ2 * pri_var) / (σ2 + nn*pri_var)
	post_mu := (pri_mu/pri_var)/(nn/σ2+1/pri_var) + samp_mean*(nn/σ2)/(nn/σ2+1/pri_var)
	post_sigma := math.Sqrt(post_var)
	return Normal_Qtl_For(post_mu, post_sigma, p)
}

// Credible interval for unknown Normal mean, with KNOWN sigma, and normal prior
// Bolstad 2007 (2e): 212, eq. 11.7
func NormMu_CrI_NPriKnown(n int, samp_mean, sigma, pri_mu, pri_sigma, alpha float64) (lo, hi float64) {
	// samp_mean		sample mean of observations taken from Normal distribution
	// sigma		standard deviation of population, assumed to be known
	// n			number of observations
	// pri_mu		Normal prior mean
	// pri_sigma		Normal prior standard deviation
	// alpha		posterior probability that the true mean lies outside the credible interval
	nn := float64(n)
	σ2 := sigma * sigma
	pri_var := pri_sigma * pri_sigma
	post_var := (σ2 * pri_var) / (σ2 + nn*pri_var)
	post_mu := (pri_mu/pri_var)/(nn/σ2+1/pri_var) + samp_mean*(nn/σ2)/(nn/σ2+1/pri_var)
//	post_mu := (pri_mu/pri_var)/(nn*samp_mean/σ2+1/pri_var) + ((nn / σ2) / (nn/σ2 + 1/pri_var))
	post_sigma := math.Sqrt(post_var)
	lo = Normal_Qtl_For(post_mu, post_sigma, alpha/2)
	hi = Normal_Qtl_For(post_mu, post_sigma, 1-alpha/2)
	return lo, hi
}

/* waiting for StudentsT_Qtl_For() to be implemented
// Credible interval for unknown Normal mean, with UNKNOWN sigma, and normal prior, equal tail area
// Bolstad 2007 (2e): 212, eq. 11.8
// n			number of observations
// samp_mean		sample mean of observations taken from Normal distribution
// samp_sigma	standard deviation of the sample
// pri_mu		Normal prior mean
// pri_sigma		Normal prior standard deviation
// alpha		posterior probability that the true mean lies outside the credible interval
// untested ...
func NormMu_CrI_NPriUnkn(n int, samp_mean, samp_sigma, pri_mu, pri_sigma, alpha float64) (lo, hi float64) {
	nn := float64(n)
	nu := float64(n - 1)
	samp_var := samp_sigma * samp_sigma
	pri_var := pri_sigma * pri_sigma
	post_var := (samp_var * pri_var) / (samp_var + nn*pri_var)
	post_mu := (pri_mu/pri_var)/(nn*samp_mean/samp_var+1/pri_var) + ((nn / samp_var) / (nn/samp_var + 1/pri_var))
	post_sigma := math.Sqrt(post_var)
	t := StudentsT_Qtl_For(alpha/2, nu)
	lo = post_mu - t*post_sigma
	hi = post_mu + t*post_sigma
	return lo, hi
}
*/

// Credible interval for unknown Normal mean, with KNOWN sigma, and flat prior
// Bolstad 2007 (2e): 212, eq. 11.7
// untested ...
// samp_mean		sample mean of observations taken from Normal distribution
// sigma		standard deviation of population, assumed to be known
// n			number of observations
// alpha		posterior probability that the true mean lies outside the credible interval
func NormMu_CrI_FPriKnown(n int, samp_mean, sigma, alpha float64) (lo, hi float64) {
	nn := float64(n)
	post_mu := samp_mean
	post_var := (sigma * sigma / nn)
	post_sigma := math.Sqrt(post_var)
	lo = Normal_Qtl_For(post_mu, post_sigma, alpha/2)
	hi = Normal_Qtl_For(post_mu, post_sigma, 1-alpha/2)
	return lo, hi
}

/* waiting for StudentsT_Qtl_For() to be implemented
// Credible interval for unknown Normal mean, with UNKNOWN sigma, and flat prior
// Bolstad 2007 (2e): 212, eq. 11.8
// samp_mean		sample mean of observations taken from Normal distribution
// sigma		standard deviation of population, unknown
// n			number of observations
// alpha		posterior probability that the true mean lies outside the credible interval
// untested ...
func NormMu_CrI_FPriUnkn(n int, samp_mean, sigma, alpha float64) (lo, hi float64) {
	nn := float64(n)
	nu := float64(n - 1)
	post_mu := samp_mean
	post_var := (sigma * sigma / nn)
	post_sigma := math.Sqrt(post_var)
	t := StudentsT_Qtl_For(alpha/2, nu)
	lo = post_mu - t*post_sigma
	hi = post_mu + t*post_sigma
	return lo, hi
}
*/

