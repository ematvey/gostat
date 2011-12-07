package stat

import (
//	"math"
	"math/rand"
)

func Uniform_PDF() func(x float64) float64 {
	return func(x float64) float64 {
		if 0 <= x && x <= 1 {
			return 1
		}
		return 0
	}
}
func Uniform_LnPDF() func(x float64) float64 {
	return func(x float64) float64 {
		if 0 <= x && x <= 1 {
			return 0
		}
		return negInf
	}
}

var NextUniform func() float64 = rand.Float64

func Uniform() func() float64 { return NextUniform }

func Exp_PDF(λ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		return λ * NextExp(-1*λ*x)
	}
}
func Exp_LnPDF(λ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return negInf
		}
		return log(λ) - λ*x
	}
}
func NextExp(λ float64) float64    { return rand.ExpFloat64() / λ }
func Exp(λ float64) func() float64 { return func() float64 { return NextExp(λ) } }

func Gamma_PDF(α float64, λ float64) func(x float64) float64 {
	expPart := Exp_PDF(λ)
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		return expPart(x) * pow(λ*x, α-1) / Γ(α)
	}
}
func Gamma_LnPDF(α float64, λ float64) func(x float64) float64 {
	expPart := Exp_LnPDF(λ)
	return func(x float64) float64 {
		if x < 0 {
			return negInf
		}
		return expPart(x) + (α-1)*log(λ*x) - LnΓ(α)
	}
}
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

func Xsquare_PDF(n int64) func(x float64) float64 {
	k := float64(n) / 2
	normalization := pow(0.5, k) / Γ(k)
	return func(x float64) float64 {
		return normalization * pow(x, k-1) * NextExp(-x/2)
	}
}
func Xsquare_LnPDF(n int64) func(x float64) float64 {
	k := float64(n) / 2
	normalization := log(0.5)*k - LnΓ(k)
	return func(x float64) float64 {
		return normalization + log(x)*(k-1) - x/2
	}
}
//Xsquare(n) => sum of n N(0,1)^2
func NextXsquare(n int64) (x float64) {
	for i := iZero; i < n; i++ {
		n := NextNormal(0, 1)
		x += n * n
	}
	return
}
func Xsquare(n int64) func() float64 {
	return func() float64 {
		return NextXsquare(n)
	}
}

func StudentsT_PDF(ν float64) func(x float64) float64 {
	normalization := Γ((ν+1)/2) / (sqrt(ν*π) * Γ(ν/2))
	return func(x float64) float64 {
		return normalization * pow(1+x*x/ν, -(ν+1)/2)
	}
}
func StudentsT_LnPDF(ν float64) func(x float64) float64 {
	normalization := LnΓ((ν+1)/2) - log(sqrt(ν*π)) - LnΓ(ν/2)
	return func(x float64) float64 {
		return normalization + log(1+x*x/ν)*-(ν+1)/2
	}
}
//StudentsT(ν) => N(0, 1)*sqrt(ν/NextGamma(ν/2, 2))
func NextStudentsT(ν float64) float64 {
	return NextNormal(0, 1) * sqrt(ν/NextGamma(ν/2, 2))
}
func StudentsT(ν float64) func() float64 {
	return func() float64 {
		return NextStudentsT(ν)
	}
}


