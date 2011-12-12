// Gamma distribution

package stat
import (
	"fmt"
	"math"
	"go-fn.googlecode.com/hg/fn"
)

// Probability density function
func Gamma_PDF(α float64, λ float64) func(x float64) float64 {
	expPart := Exp_PDF(λ)
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		return expPart(x) * pow(λ*x, α-1) / Γ(α)
	}
}

// Another algorithm for Gamma PDF
func Gamma_PDF_2(k float64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		return pow(x, k-1) * exp(-x/θ) / (Γ(k) * pow(θ, k))
	}
}

// Natural logarithm of the probability density function
func Gamma_LnPDF(α float64, λ float64) func(x float64) float64 {
	expPart := Exp_LnPDF(λ)
	return func(x float64) float64 {
		if x < 0 {
			return negInf
		}
		return expPart(x) + (α-1)*log(λ*x) - LnΓ(α)
	}
}

// Random value drawn from the distribution
func NextGamma(α float64, λ float64) float64 {
	//if α is a small integer, this way is faster on my laptop
	if α == float64(int64(α)) && α <= 15 {
		x := NextExp(λ)
		for i := 1; i < int(α); i++ {
			x += NextExp(λ)
		}
		return x
	}

	if α < 0.75 {
		return RejectionSample(Gamma_PDF(α, λ), Exp_PDF(λ), Exp(λ), 1)
	}

	//Tadikamalla ACM '73
	A := α - 1
	B := 0.5 + 0.5*sqrt(4*α-3)
	C := A * (1 + B) / B
	D := (B - 1) / (A * B)
	S := A / B
	P := 1.0 / (2 - exp(-S))
	var x, y float64
	for i := 1; ; i++ {
		u := NextUniform()
		if u > P {
			var E float64
			for E = -log((1 - u) / (1 - P)); E > S; E = E - A/B {
			}
			x = A - B*E
			y = A - x
		} else {
			x = A - B*log(u/P)
			y = x - A
		}
		u2 := NextUniform()
		if log(u2) <= A*log(D*x)-x+y/B+C {
			break
		}
	}
	return x / λ
}

func Gamma(α float64, λ float64) func() float64 {
	return func() float64 { return NextGamma(α, λ) }
}

// Cumulative distribution function
func Gamma_CDF(k float64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if k < 0 || θ < 0 {
			panic(fmt.Sprintf("k < 0 || θ < 0"))
		}
		if x < 0 {
			return 0
		}
		return fn.Iγ(k, x/θ) / fn.Γ(k)
	}
}

// Value of the probability density function at x
func Gamma_PDF_At(k , θ, x float64)  float64 {
	pdf := Gamma_PDF(k , θ)
	return pdf(x)
}

// Value of the cumulative distribution function at x
func Gamma_CDF_At(k , θ, x float64)  float64 {
	cdf := Gamma_CDF(k , θ)
	return cdf(x)
}

// Inverse CDF (Quantile) function
func Gamma_InvCDF(k float64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		var eps, y_new, h float64
		eps = 1e-4
		y := k * θ
		y_old := y
	L:
		for i := 0; i < 100; i++ {
			h = (Gamma_CDF_At(k, θ, y_old) - x) / Gamma_PDF_At(k, θ, y_old)
			y_new = y_old - h
			if y_new <= eps {
				y_new = y_old / 10
				h = y_old - y_new
			}
			if math.Abs(h) < eps {
				break L
			}
			y_old = y_new
		}
		return y_new
	}
}

// Value of the inverse CDF for probability p
func Gamma_InvCDF_For(k, θ, p float64)  float64 {
	cdf:=Gamma_InvCDF(k, θ)
	return cdf(p)
}


