// Bayesian inference about the params (μ, σ) of Normal (Gaussian) distribution.
// Bolstad 2007 (2e): Chapter 11, p. 199 and further.

package bayes

import (
	"fmt"
	. "gostat.googlecode.com/hg/stat/prob"
	"math"
)

// PMF of the posterior distribution of unknown Normal μ, with KNOWN σ, and discrete prior, for single observation. 
// Bolstad 2007 (2e): 200-201.
func NormMuSingle_PMF_DPri(y, σ float64, μ []float64, μPri []float64) (post []float64) {
// y	single observation taken from Normal distribution
// σ	standard deviation of population, assumed to be known
// μ	array of possible discrete values of μ
// μPri	array of associated prior probability masses
	nPoss := len(μ)
	if len(μPri) != nPoss {
		panic(fmt.Sprintf("len(μ) != len(μPri)"))
	}
	post = make([]float64, nPoss)
	sum := 0.0
	for i := 0; i< nPoss; i++ {
		z := (y-μ[i])/σ
		like := Z_PDF_At(z)
		post[i] = μPri[i] * like
		sum += post[i]
	}
	for i := 0; i< nPoss; i++ {
		post[i] /= sum
	}
	return
}

// PMF of the posterior distribution of unknown Normal μ, with KNOWN σ, and discrete prior, for sample
// Bolstad 2007 (2e): 203, eq. 11.2
func NormMu_PMF_DPri(nObs int, ȳ, σ float64, μ []float64, μPri []float64) (post []float64) {
// nObs		number of observations in the sample (= length of the samle array)
// ȳ		sample mean of the observed values
// σ		standard deviation of population, assumed to be known
// μ		array of possible discrete values of μ
// μPri		array of associated prior probability masses
	nPoss := len(μ)	// number of possible values of the parameter μ
	if len(μPri) != nPoss {
		panic(fmt.Sprintf("len(μ) != len(μPri)"))
	}
	post = make([]float64, nPoss)
	n := float64(nObs)
	sum := 0.0
	for i := 0; i< nPoss; i++ {
		σ2 := σ * σ
		ẟ:=ȳ-μ[i]
		like := math.Exp(-1/(2*σ2/n)*ẟ*ẟ)
		post[i] = μPri[i] * like
		sum += post[i]
	}
	for i := 0; i< nPoss; i++ {
		post[i] /= sum
	}
	return
}

// Posterior mean for unknown Normal μ, with KNOWN σ. 
// Bolstad 2007 (2e): 209, eq. 11.6
func NormMuPostMean(nObs int, ȳ, σ, μPri, σPri float64) float64 {
	// μPri		prior mean
	// σPri		prior standard deviation
	// nObs			size of sample == number of measurements
	// σ		standard deviation of population, assumed to be known (alternatively, use an estimate)
	σ2 := σ * σ
	n := float64(nObs)
	σ2Pri := σPri * σPri
	μPost := (μPri/σ2Pri)/(n/σ2+1/σ2Pri) + ȳ*(n/σ2)/(n/σ2+1/σ2Pri)
	return (μPost)
}

// Posterior standard deviation for unknown Normal μ, with KNOWN σ. 
// Bolstad 2007 (2e): 209, eq. 11.5
func NormMuPostStd(nObs int, σ, μPri, σPri float64) float64 {
	// μPri		prior mean
	// σPri		prior standard deviation
	// nObs		size of sample == number of measurements
	// σ		standard deviation of population, assumed to be known (alternatively, use an estimate)
	σ2 := σ * σ
	n := float64(nObs)
	σ2Pri := σPri * σPri
	σ2Post := (σ2 * σ2Pri) / (σ2 + n*σ2Pri)
	σPost := math.Sqrt(σ2Post)
	return (σPost)
}

// Quantile for posterior distribution of unknown Normal μ, with KNOWN σ, and flat prior (Jeffrey's prior), for single observation
// Bolstad 2007 (2e): 206
func NormMuSingle_Qtl_FPri(y, σ, p float64) float64 {
// y		single observation taken from Normal distribution
// σ		standard deviation of population, assumed to be known
// p		probability for which the quantile will be returned
// untested ...
	μPost := y
	σPost := σ
	qtl := Normal_Qtl_For(μPost, σPost, p)
	return (qtl)
}

// Quantile for posterior distribution of unknown Normal μ, with KNOWN σ, and flat prior (Jeffrey's prior), for sample
// Bolstad 2007 (2e): 207
func NormMu_Qtl_FPri(nObs int, ȳ, σ, p float64) float64 {
// ȳ		sample mean of observations taken from Normal distribution
// σ		standard deviation of population, assumed to be known
// nObs		number of observations
// p		probability for which the quantile will be returned

	if σ <= 0 {
		panic(fmt.Sprintf("Prior standard deviation must be greater than zero"))
	}

	n := float64(nObs)
	σ2 := σ * σ
	μPost := ȳ
	σ2Post := σ2 / n
	σPost := math.Sqrt(σ2Post)
	return Normal_Qtl_For(μPost, σPost, p)
}

// Quantile for posterior distribution of unknown Normal μ, with KNOWN σ, and Normal prior, for single observation
// Bolstad 2007 (2e): 208, eq. 11.4
func NormMuSingle_Qtl_NPri(y, σ, μPri, σPri, p float64) float64 {
// y		single observation taken from Normal distribution
// σ	standard deviation of population, assumed to be known
// μPri	Normal prior mean
// σPri	Normal prior standard deviation
// p		probability for which the quantile will be returned
// untested ...
	if σ <= 0 {
		panic(fmt.Sprintf("Prior standard deviation must be greater than zero"))
	}

	σ2 := σ * σ
	σ2Pri := σPri * σPri
	μPost := (σ2*μPri + σ2Pri*y) / (σ2 + σ2Pri)
	σ2Post := (σ2 * σ2Pri) / (σ2 + σ2Pri)
	σPost := math.Sqrt(σ2Post)
	return Normal_Qtl_For(μPost, σPost, p)
}

// Quantile for posterior distribution of unknown Normal μ, with KNOWN σ, and Normal prior, for sample
// Bolstad 2007 (2e): 209, eq. 11.5, 11.6
func NormMu_Qtl_NPri(nObs int, ȳ, σ, μPri, σPri, p float64) float64 {
// ȳ		sample mean of observations taken from Normal dist
// σ		standard deviation of population, assumed to be known
// nObs			number of observations
// μPri		Normal prior mean
// σPri		Normal prior standard deviation
// p			probability for which the quantile will be returned
	n := float64(nObs)
	σ2 := σ * σ
	σ2Pri := σPri * σPri
	σ2Post := (σ2 * σ2Pri) / (σ2 + n*σ2Pri)
	μPost := (μPri/σ2Pri)/(n/σ2+1/σ2Pri) + ȳ*(n/σ2)/(n/σ2+1/σ2Pri)
	σPost := math.Sqrt(σ2Post)
	return Normal_Qtl_For(μPost, σPost, p)
}

// Credible interval for unknown Normal μ, with KNOWN σ, and Normal prior
// Bolstad 2007 (2e): 212, eq. 11.7
func NormMu_CrI_NPriKnown(nObs int, ȳ, σ, μPri, σPri, α float64) (lo, hi float64) {
	// ȳ		sample mean of observations taken from Normal distribution
	// σ		standard deviation of population, assumed to be known
	// nObs			number of observations
	// μPri		Normal prior mean
	// σPri		Normal prior standard deviation
	// α		posterior probability that the true μ lies outside the credible interval
	n := float64(nObs)
	σ2 := σ * σ
	σ2Pri := σPri * σPri
	σ2Post := (σ2 * σ2Pri) / (σ2 + n*σ2Pri)
	μPost := (μPri/σ2Pri)/(n/σ2+1/σ2Pri) + ȳ*(n/σ2)/(n/σ2+1/σ2Pri)
//	μPost := (μPri/σ2Pri)/(n*ȳ/σ2+1/σ2Pri) + ((n / σ2) / (n/σ2 + 1/σ2Pri))
	σPost := math.Sqrt(σ2Post)
	lo = Normal_Qtl_For(μPost, σPost, α/2)
	hi = Normal_Qtl_For(μPost, σPost, 1-α/2)
	return lo, hi
}

/* waiting for StudentsT_Qtl_For() to be implemented
// Credible interval for unknown Normal μ, with UNKNOWN σ, and Normal prior, equal tail area
// Bolstad 2007 (2e): 212, eq. 11.8
func NormMu_CrI_NPriUnkn(nObs int, ȳ, samp_σ, μPri, σPri, α float64) (lo, hi float64) {
// nObs			number of observations
// ȳ		sample mean of observations taken from Normal distribution
// samp_σ	standard deviation of the sample
// μPri		Normal prior mean
// σPri		Normal prior standard deviation
// α		posterior probability that the true μ lies outside the credible interval
// untested ...
	n := float64(nObs)
	nu := float64(nObs - 1)
	samp_var := samp_σ * samp_σ
	σ2Pri := σPri * σPri
	σ2Post := (samp_var * σ2Pri) / (samp_var + n*σ2Pri)
	μPost := (μPri/σ2Pri)/(n*ȳ/samp_var+1/σ2Pri) + ((n / samp_var) / (n/samp_var + 1/σ2Pri))
	σPost := math.Sqrt(σ2Post)
	t := StudentsT_Qtl_For(α/2, nu)
	lo = μPost - t*σPost
	hi = μPost + t*σPost
	return lo, hi
}
*/

// Credible interval for unknown Normal μ, with KNOWN σ, and flat prior
// Bolstad 2007 (2e): 212, eq. 11.7
func NormMu_CrI_FPriKnown(nObs int, ȳ, σ, α float64) (lo, hi float64) {
// untested ...
// ȳ		sample mean of observations taken from Normal distribution
// σ		standard deviation of population, assumed to be known
// nObs		number of observations
// α		posterior probability that the true μ lies outside the credible interval
	n := float64(nObs)
	μPost := ȳ
	σ2Post := (σ * σ / n)
	σPost := math.Sqrt(σ2Post)
	lo = Normal_Qtl_For(μPost, σPost, α/2)
	hi = Normal_Qtl_For(μPost, σPost, 1-α/2)
	return lo, hi
}

/* waiting for StudentsT_Qtl_For() to be implemented
// Credible interval for unknown Normal μ, with UNKNOWN σ, and flat prior
// Bolstad 2007 (2e): 212, eq. 11.8
func NormMu_CrI_FPriUnkn(nObs int, ȳ, σ, α float64) (lo, hi float64) {
// ȳ		sample mean of observations taken from Normal distribution
// σ		standard deviation of population, unknown
// nObs		number of observations
// α		posterior probability that the true μ lies outside the credible interval
// untested ...
	n := float64(nObs)
	nu := float64(nObs - 1)
	μPost := ȳ
	σ2Post := (σ * σ / n)
	σPost := math.Sqrt(σ2Post)
	t := StudentsT_Qtl_For(α/2, nu)
	lo = μPost - t*σPost
	hi = μPost + t*σPost
	return lo, hi
}
*/

