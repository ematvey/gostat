// Pareto Type II Distribution
// params:
// θ > 0.0 (scale)
// α > 0.0 (shape)
// support: x >= θ
// inspired by R:actuar

package prob

import (
	"math"
	"math/rand"
	. "go-fn.googlecode.com/hg/fn"
)

func ParetoII_ChkParams(θ, α float64) bool {
	ok := true
	if α <= 0 || θ <= 0  {
		ok = false
	}
	return ok
}

func ParetoII_ChkSupport(x float64) bool {
	ok := true
	if x < 0 {
		ok = false
	}
	return ok
}

func ParetoII_PDF(θ, α float64) func(x float64) float64 {
	// We work with the density expressed as
	// α * u^α * (1 - u) / x
	// with u = 1/(1 + v), v = x/θ.
	return func(x float64) float64 {
		var p float64
		if x < 0 {
			p = 0
		} else if x == 0 {
			p = α / θ
		} else {
		tmp := math.Log(x) - math.Log(θ)
		logu := - math.Log1p(math.Exp(tmp))
		log1mu := - math.Log1p(math.Exp(-tmp))
		p = math.Exp(math.Log(α) + α * logu + log1mu - math.Log(x))
		}
		return p
	}
}

func ParetoII_PDF_At(θ, α, x float64) float64 {
	pdf := ParetoII_PDF(θ, α)
	return pdf(x)
}

func ParetoII_CDF(θ, α float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		u := math.Exp(-math.Log1p(math.Exp(math.Log(x) - math.Log(θ))))
		return 1 - math.Pow(u, α)
	}
}

func ParetoII_CDF_At(θ, α, q, x float64) float64 {
	cdf := ParetoII_CDF(θ, α)
	return cdf(x)
}


// Inverse of the cumulative Pareto Type II probability density function (quantile).
func ParetoII_Qtl(θ, α float64) func(p float64) float64 {
	return func(p float64) float64 {
		if p < 0 || p > 1 {
			panic("probability out of range 0..1")
		}
		return θ * (math.Pow((0.5 - (p) + 0.5), -1.0 / α) - 1.0)
	}
}

// Inverse of the cumulative Pareto Type II probability density function (quantile) for given probability.
func ParetoII_Qtl_For(θ, α, p float64) float64 {
	cdf := Pareto_Qtl(θ, α)
	return cdf(p)
}

func NextParetoII(θ, α float64) float64 {
	qtl := ParetoII_Qtl(θ, α)
	p := rand.Float64()
	return qtl(p)
}

func ParetoII_Moment(θ, α float64, order int) float64 {
	o := float64(order)
	if α <= o {
		panic("not defined")
	}
	return math.Pow(θ, o) * Γ(1.0 + o) * Γ(α - o) / Γ(α)
}

func ParetoII_Mean(θ, α float64) float64 {
	return ParetoII_Moment(θ, α, 1)
}

func ParetoII_Var(θ, α float64) float64 {
	return ParetoII_Moment(θ, α, 2)
}

func ParetoII_Skew(θ, α float64) float64 {
	return ParetoII_Moment(θ, α, 3)
}

