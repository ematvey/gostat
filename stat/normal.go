package stat

import (
	"math"
	"math/rand"
)

func Normal_PDF(μ float64, σ float64) func(x float64) float64 {
	normal_normalizer := 0.3989422804014327 / σ
	return func(x float64) float64 { return normal_normalizer * exp(-1*(x-μ)*(x-μ)/(2*σ*σ)) }
}
func Normal_LnPDF(μ float64, σ float64) func(x float64) float64 {
	ln_normal_normalizer := -0.91893853320467267 - log(σ)
	return func(x float64) float64 { return ln_normal_normalizer - (x-μ)*(x-μ)/(2*σ*σ) }
}
func NextNormal(μ float64, σ float64) float64 { return rand.NormFloat64()*σ + μ }
func Normal(μ , σ float64) func() float64 {
	return func() float64 { return NextNormal(μ, σ) }
}

// Cumulative Distribution Function for the Normal distribution
func Normal_CDF(μ, σ float64) func(x float64) float64 {
	return func(x float64) float64 { return ((1.0 / 2.0) * (1 + math.Erf((x-μ)/(σ*math.Sqrt2)))) }
}

// Inverse CDF of standardized Normal distribution for probability p
func ZInv_CDF_For(p float64) float64 {
	return 0	// to be implemented
}


