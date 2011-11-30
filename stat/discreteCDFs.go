package stat
import "math"
func Binomial_CDF_trivial(ρ float64, n int64) func(k int64) float64 {  // trivial (but working) implementation, redo with Incomplete Beta
	return func(k int64) float64 {
		var p float64 = 0
		var i int64
		pmf:=Binomial_PMF(ρ , n)
			for i = 0; i<=k; i++ {
				p+=pmf(i)
			}
		return p
	}
}

func Binomial_CDF(ρ float64, n int64) func(k int64) float64 {
	return func(k int64) float64 {
		p:= Beta_CDF_At((float64)(n-k), (float64)(k+1), 1-ρ)
		return p
	}
}

func Poisson_CDF_trivial(λ float64) func(k int64) float64 {  // trivial (not working due to bug in Poisson_PMF(λ) implementation, redo with Incomplete Gamma
	return func(k int64) float64 {
		var p float64 = 0
		var i int64
		pmf:=Poisson_PMF(λ)
			for i = 0; i<=k; i++ {
				p+=pmf(i)
			}
		return p
	}
}
func Poisson_CDF(λ float64) func(k int64) float64 { 
	return func(k int64) float64 {
		p:=math.Exp(math.Log(IΓ(((float64)(k+1)), λ)) - (LnFact(k)))
		return p
	}
}

func LnPoisson_CDF(λ float64) func(k int64) float64 { 
	return func(k int64) float64 {
		k1:=(float64)(k+1)
		return log(IΓ(k1, λ)) - LnFact(k)
	}
}


func NegativeBinomial_CDF(ρ float64, r int64) func(k int64) float64 {
	return func(k int64) float64 {
		Ip:=Beta_CDF_At(float64(k+1), float64(r), ρ)
		return 1-Ip
	}
}

