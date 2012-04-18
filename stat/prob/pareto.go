// Pareto Type I Distribution
// params:
// θ > 0.0 (scale)
// α > 0.0 (shape)
// support: x >= θ

package prob

import (
	"math"
)

func Pareto_PDF(θ, α float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < θ {
			return 0
		}
		t1 := α*math.Pow(θ, α)
		t2 := math.Pow(x, α + 1)
		return t1/t2
	}
}

func Pareto_PDF_At(θ, α, x float64) float64 {
	pdf := Pareto_PDF(θ, α)
	return pdf(x)
}

func Pareto_CDF(θ, α float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < θ {
			return 0
		}
		return 1-math.Pow(θ/x, α)
	}
}

func Pareto_CDF_At(θ, α, x float64) float64 {
	cdf := Pareto_CDF(θ, α)
	return cdf(x)
}


func Pareto_Mean(θ, α float64) float64 {
	if α <= 1 {
		panic("not defined")
	}
	return α*θ/(α-1)
}

func Pareto_Median(θ, α float64) float64 {
	return θ * math.Pow(2, 1/α)
}

func Pareto_Mode(θ, α float64) float64 {
	return θ
}

func Pareto_Var(θ, α float64) float64 {
	if α <= 2 {
		panic("not defined")
	}
	α1:= (α-1)
	α2:= (α-2)
	return θ*θ*α/(α1*α1*α2)
}


