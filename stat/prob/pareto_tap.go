// Tapered Pareto Distribution
// params: 
// θ > 0.0 (scale) "a"
// α > 0.0 (shape) "lambda"
// taper > 0.0 (tapering) 
// support: x >= θ 
// inspired by R:PtProcess

package prob

import (
	"math"
	"math/rand"
)
func ParetoTap_ChkParams(θ, α, taper float64) bool {
	ok := true
	if α <= 0 || θ <= 0 || taper <= 0 {
		ok = false
	}
	return ok
}

func ParetoTap_ChkSupport(x float64) bool {
	ok := true
	if x < 0 {
		ok = false
	}
	return ok
}

func ParetoTap_PDF(θ, α, taper float64) func(x float64) float64 {
	return func(x float64) float64 {
		return (α/x + 1/taper) * math.Pow((θ/x), α) * math.Exp((θ-x)/taper)
	}
}

func ParetoTap_PDF_At(θ, α, taper, x float64) float64 {
	pdf := ParetoTap_PDF(θ, α, taper)
	return pdf(x)
}

func ParetoTap_CDF(θ, α, taper float64) func(x float64) float64 {
	return func(x float64) float64 {
		return 1 - math.Pow(θ/x, α) * math.Exp((θ-x)/taper)
	}
}

func ParetoTap_CDF_At(θ, α, taper, x float64) float64 {
	cdf := ParetoTap_CDF(θ, α, taper)
	return cdf(x)
}

// Inverse of the cumulative Pareto Type I probability density function (quantile).
func ParetoTap_Qtl(θ, α, taper float64) func(p float64) float64 {
	return func(p float64) float64 {
		tol:=1e-8
		x := θ+1

		// solve using Newton-Raphson method
		for {
        		delta := (ParetoTap_CDF_At(θ, α, taper, x) - p)/ ParetoTap_PDF_At(θ, α, taper, x)
			x -= delta
        		if math.Abs(delta) < tol {
				break
			}
		}
		return x
	}
}

// Inverse of the cumulative Pareto Type I probability density function (quantile) for given probability.
func ParetoTap_Qtl_For(θ, α, taper, p float64) float64 {
	cdf := ParetoTap_Qtl(θ, α, taper)
	return cdf(p)
}

func NextParetoTap(θ, α, taper float64) float64 {
	qtl := ParetoTap_Qtl(θ, α, taper)
	p := rand.Float64()
	return qtl(p)
}

/*

func ParetoTap_Mean(θ, α, taper float64) float64 {
}

func ParetoTap_Median(θ, α, taper float64) float64 {
}

func ParetoTap_Mode(θ, α, taper float64) float64 {
}

func ParetoTap_Var(θ, α, taper float64) float64 {
}

*/



