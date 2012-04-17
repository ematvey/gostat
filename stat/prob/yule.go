// Yule–Simon distribution
// Yule, G. U. (1925). "A Mathematical Theory of Evolution, based on the Conclusions of Dr. J. C. Willis, F.R.S". Philosophical Transactions of the Royal Society of London, Ser. B 213 (402–410): 21–87. doi:10.1098/rstb.1925.0002
// Simon, H. A. (1955). "On a class of skew distribution functions". Biometrika 42 (3–4): 425–440. doi:10.1093/biomet/42.3-4.425

package prob

import (
	"math"
	. "go-fn.googlecode.com/hg/fn"
)

// Probability Mass Function for the Yule distribution
func YulePMF(a float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := a*B(a+1, float64(k))
		return p
	}
}

// Probability Mass Function for the Yule distribution at k
func YulePMF_At(a float64, k int64) float64 {
	pmf := YulePMF(a)
	return pmf(k)
}

// Cumulative Distribution Function for the Yule distribution
func YuleCDF(a float64) func(k int64) float64 {
	return func(k int64) float64 {
		kk := float64(k)
		p := 1- kk*B(kk, a+1)
		return p
	}
}

// Cumulative Distribution Function for the Yule distribution at k
func YuleCDF_At(a float64, k int64) float64 {
	cdf := YuleCDF(a)
	return cdf(k)
}

// Yule distributed random variate
// Devroye 1986: 553.
// Devroye, L. 1986: Non-Uniform Random Variate Generation. Springer-Verlag, New York. ISBN 0-387-96305-7.

func NextYule(a float64) (k int64) {
	e1 := NextExp(2)
	e2 := NextExp(2)
	k = int64(math.Ceil(-e1 / (math.Log(1-math.Exp(-e2/(a-1))))))
        return 
}

func Yule(a float64) func() int64 {
	return func() int64 { return NextYule(a) }
}

func YuleMean(a float64) float64 {
        if a <=1 {
		panic("not defined")
	}
        return a/(a-1)
}

func YuleMode() float64 {
        return 1.00
}

func YuleVar(a float64) float64 {
        if a <=2 {
		panic("not defined")
	}
	aa := a*a
	a1 := a-1
	a2 := a-2
        return aa/(a1*a1*a2)
}

