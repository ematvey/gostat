package stat

import (
	"math"
	. "go-fn.googlecode.com/hg/fn"
)

// Probability Mass Function for the Hypergeometric distribution
func Hypergeometric_PMF(size, m, n int64) func(i int64) float64 {
	return func(k int64) float64 {
		if size < 1 || m < 0 || m > size || n < 1 || n > size {
			panic("bad param")
		} 
		p := BinomCoeff(m, k) * BinomCoeff(size-m, n-k)  / BinomCoeff(size, n) 
		return p
	}
}

func Hypergeometric_PMF_At(size, m, n, k int64) float64 {
	if float64(k) < math.Max(0, float64(n+m-size)) || float64(k) > math.Min(float64(m), float64(n)) {
		panic("bad param")
	} 
	pmf := Hypergeometric_PMF(size, m, n)
	return pmf(k)
}

// Cumulative Distribution Function for the Hypergeometric distribution
func Hypergeometric_CDF(size, m, n int64) func(k int64) float64 {
	return func(k int64) float64 {
		var p float64 = 0
		var i int64
		pmf:=Hypergeometric_PMF(size, m, n)
			for i = 0; i<=k; i++ {
				p+=pmf(i)
			}
		return p
	}
}

func Hypergeometric_CDF_At(size, m, n, k int64) float64 {
	cdf := Hypergeometric_CDF(size, m, n)
	return cdf(k)
}


