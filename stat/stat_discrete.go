package stat

import (
	"math"
	"math/rand"
)

func Range_PMF(n int64) func(i int64) float64 {
	return func(i int64) float64 {
		return fOne / float64(n)
	}
}
func LnRange_PMF(n int64) func(i int64) float64 {
	return func(i int64) float64 {
		return -log(float64(n))
	}
}
func NextRange(n int64) int64 {
	return rand.Int63n(n)
}
func Range(n int64) func() int64 {
	return func() int64 {
		return NextRange(n)
	}
}

func Choice_PMF(θ []float64) func(i int64) float64 {
	return func(i int64) float64 {
		return θ[i]
	}
}
func Choice_LnPMF(θ []float64) func(i int64) float64 {
	return func(i int64) float64 {
		return log(θ[i])
	}
}
func NextChoice(θ []float64) int64 {
	u := NextUniform()
	i := 0
	sum := θ[0]
	for ; sum < u && i < len(θ)-1; i++ {
		sum += θ[i+1]
	}
	if u >= sum {
		return int64(len(θ))
	}
	return int64(i)
}
func Choice(θ []float64) func() int64 {
	return func() int64 {
		return NextChoice(θ)
	}
}
func NextLogChoice(lws []float64) int64 {
	return LogChoice(lws)()
}
func LogChoice(lws []float64) func() int64 {
	max := lws[0]
	for _, lw := range lws[1:len(lws)] {
		if lw > max {
			max = lw
		}
	}
	ws := make([]float64, len(lws))
	var sum float64
	for i, lw := range lws {
		ws[i] = math.Exp(lw - max)
		sum += ws[i]
	}
	norm := 1 / sum
	for i := range ws {
		ws[i] *= norm
	}
	return Choice(ws)
}
func Multinomial_PMF(θ []float64, n int64) func(x []int64) float64 {
	return func(x []int64) float64 {
		if len(x) != len(θ) {
			return 0
		}
		l := fOne
		totalx := iZero
		for i := 0; i < len(x); i++ {
			l *= pow(θ[i], float64(x[i]))
			l /= Γ(float64(x[i] + 1))
			totalx += x[i]
		}
		if totalx != n {
			return 0
		}
		l *= Γ(float64(totalx + 1))
		return l
	}
}
func Multinomial_LnPMF(θ []float64, n int64) func(x []int64) float64 {
	return func(x []int64) float64 {
		if len(x) != len(θ) {
			return negInf
		}
		l := fZero
		totalx := iZero
		for i := 0; i < len(x); i++ {
			l += log(θ[i]) * float64(x[i])
			l -= LnΓ(float64(x[i] + 1))
			totalx += x[i]
		}
		if totalx != n {
			return negInf
		}
		l += LnΓ(float64(totalx + 1))
		return l
	}
}
func NextMultinomial(θ []float64, n int64) []int64 {
	x := make([]int64, len(θ))
	chooser := Choice(θ)
	for i := iZero; i < n; i++ {
		x[chooser()]++
	}
	return x
}
func Multinomial(θ []float64, n int64) func() []int64 {
	return func() []int64 {
		return NextMultinomial(θ, n)
	}
}

func Bernoulli_PMF(ρ float64) func(i int64) float64 {
	return func(x int64) float64 {
		if x == 1 {
			return ρ
		}
		return 1 - ρ
	}
}
func Bernoulli_LnPMF(ρ float64) func(i int64) float64 {
	return func(x int64) float64 {
		if x == 1 {
			return log(ρ)
		}
		return log(1 - ρ)
	}
}
func NextBernoulli(ρ float64) int64 {
	if NextUniform() < ρ {
		return 1
	}
	return 0
}
func Bernoulli(ρ float64) func() int64 { return func() int64 { return NextBernoulli(ρ) } }

func Geometric_PMF(ρ float64) func(i int64) float64 {
	return func(n int64) float64 { return ρ * pow(ρ, float64(n)) }
}
func Geometric_LnPMF(ρ float64) func(i int64) float64 {
	return func(n int64) float64 { return log(ρ) + float64(n)*log(ρ) }
}
//NextGeometric(ρ) => # of NextBernoulli(ρ) failures before one success
func NextGeometric(ρ float64) int64 {
	if NextBernoulli(ρ) == 1 {
		return 1 + NextGeometric(ρ)
	}
	return 0
}
func Geometric(ρ float64) func() int64 { return func() int64 { return NextGeometric(ρ) } }

func Binomial_PMF(ρ float64, n int64) func(i int64) float64 {
	return func(i int64) float64 {
		p := pow(ρ, float64(i)) * pow(1-ρ, float64(n-i))
		p *= Γ(float64(n+1)) / (Γ(float64(i+1)) * Γ(float64(n-i+1)))
		return p
	}
}
func Binomial_LnPMF(ρ float64, n int64) func(i int64) float64 {
	return func(i int64) float64 {
		p := log(ρ)*float64(i) + log(1-ρ)*float64(n-i)
		p += LnΓ(float64(n+1)) - LnΓ(float64(i+1)) - LnΓ(float64(n-i+1))
		return p
	}
}
func NextBinomial(ρ float64, n int64) (result int64) {
	for i := int64(0); i <= n; i++ {
		result += NextBernoulli(ρ)
	}
	return
}
func Binomial(ρ float64, n int64) func() int64 {
	return func() int64 { return NextBinomial(ρ, n) }
}

/*
func Poisson_LnPMF(λ float64) (foo func(i int64) float64) {
	pmf := Poisson_PMF(λ)
	return func(i int64) (p float64) {
		return log(pmf(i))
		//p = -λ +log(λ)*float64(i)
		//x := log(Γ(float64(i)+1))
		//_ = x
		//p -= LnΓ(float64(i)+1)
		//return p
	}
}
*/
func Poisson_LnPMF(λ float64) func(k int64) float64 {
	return func(k int64) (p float64) {
		i := float64(k)
		a := log(λ) * i
		b := log(Γ(i + 1))
		p = a - b - λ
		return p
	}
}

/*
func Poisson_PMF(λ float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := NextExp(-λ) * pow(λ, float64(k)) / Γ(float64(k)+1)
		return p
	}
}

func Poisson_PMF(λ float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := math.Exp(-λ) * pow(λ, float64(k)) / Γ(float64(k)+1)
		return p
	}
}
*/

func Poisson_PMF(λ float64) func(k int64) float64 {
	pmf := Poisson_LnPMF(λ)
	return func(k int64) float64 {
		p := math.Exp(pmf(k))
		return p
	}
}

func NextPoisson(λ float64) int64 {
	// this can be improved upon
	i := iZero
	t := exp(-λ)
	p := fOne
	for ; p > t; p *= NextUniform() {
		i++
	}
	return i
}
func Poisson(λ float64) func() int64 {
	return func() int64 {
		return NextPoisson(λ)
	}
}

func NegativeBinomial_PMF(ρ float64, r int64) func(i int64) float64 {
	return func(k int64) float64 {
		return float64(Choose(k+r-1, r-1)) * pow(ρ, float64(r)) * pow(1-ρ, float64(k))
	}
}
func NegativeBinomial_LnPMF(ρ float64, r int64) func(i int64) float64 {
	return func(k int64) float64 {
		return LnChoose(k+r-1, r-1) + log(ρ)*float64(r) + log(1-ρ)*float64(k)
	}
}
//NegativeBinomial(ρ, r) => number of NextBernoulli(ρ) failures before r successes
func NextNegativeBinomial(ρ float64, r int64) int64 {
	k := iZero
	for r >= 0 {
		i := NextBernoulli(ρ)
		r -= i
		k += (1 - i)
	}
	return k
}
func NegativeBinomial(ρ float64, r int64) func() int64 {
	return func() int64 {
		return NextNegativeBinomial(ρ, r)
	}
}
