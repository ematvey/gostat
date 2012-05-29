// Pareto single-parameter Distribution
// The single-parameter Pareto distribution used in these functions has cumulative distribution function:
//
//	Pr[X <= x] = 1 - (μ/x)^shape, x > 0.
//
// See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second Edition, Wiley, 2004.
//
// params: 
// α > 0.0 (shape) 
// μ > 0.0 (minimum) 
// support: x >= θ 
// inspired by R:actuar

package prob

import (
	"math"
	"math/rand"
)

func ParetoSing_ChkParams(α, μ float64) bool {
	ok := true
	if α <= 0 || μ <= 0  {
		ok = false
	}
	return ok
}

func ParetoSing_ChkSupport(x, μ float64) bool {
	ok := true
	if x < μ {
		ok = false
	}
	return ok
}

func ParetoSing_PDF(α, μ float64) func(x float64) float64 {
	return func(x float64) float64 {
		p:= 0.0
		if x >= μ {
			p = math.Exp(math.Log(α) + α * math.Log(μ) - (α + 1.0) * math.Log(x))
		}
		return p
	}
}

func ParetoSing_PDF_At(α, μ, x float64) float64 {
	pdf := ParetoSing_PDF(α, μ)
	return pdf(x)
}

func ParetoSing_CDF(α, μ float64) func(x float64) float64 {
	return func(x float64) float64 {
		p:= 0.0
		if x > μ {
			p = (0.5 - (math.Pow(μ / x, α)) + 0.5)
		}
		return p
	}
}

func ParetoSing_CDF_At(α, μ, x float64) float64 {
	cdf := ParetoSing_CDF(α, μ)
	return cdf(x)
}

// Inverse of the cumulative single-parameter Pareto probability density function (quantile).
func ParetoSing_Qtl(α, μ float64) func(p float64) float64 {
	return func(p float64) float64 {
		return μ / math.Pow((0.5 - (p) + 0.5), 1.0 / α)
	}
}

// Inverse of the cumulative single-parameter Pareto probability density function (quantile) for given probability.
func ParetoSing_Qtl_For(α, μ, p float64) float64 {
	cdf := ParetoSing_Qtl(α, μ)
	return cdf(p)
}

func NextParetoSing(α, μ float64) float64 {
	qtl := ParetoSing_Qtl(α, μ)
	p := rand.Float64()
	return qtl(p)
}

func ParetoSing_Moment(α, μ float64, order int) float64 {
	o:=float64(order)
    if o >= α {
	return math.Inf(+1)
	}

    return α * math.Pow(μ, o) / (α - o)
}


