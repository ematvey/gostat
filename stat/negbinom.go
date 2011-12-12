// Negative Binomial distribution

package stat

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

func NegativeBinomial_CDF(ρ float64, r int64) func(k int64) float64 {
	return func(k int64) float64 {
		Ip:=Beta_CDF_At(float64(k+1), float64(r), ρ)
		return 1-Ip
	}
}

