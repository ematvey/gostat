package stat

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


