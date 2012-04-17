// Zeta distribution
// s > 1.0 && k > 0
package prob

import (
	"math"
	"math/rand"
	. "go-fn.googlecode.com/hg/fn"
)
// The generalized harmonic number of order n of m is given by
func H(n int64, m float64) float64 {
	var i int64
	h := 0.0
	for i = 1; i <= n; i++ {
		h += math.Pow(float64(i), m)
	}
	return h
}


// Probability Mass Function for the Zeta distribution
func Zeta_PMF(s float64) func(k int64) float64 {
	return func(k int64) float64 {
		t1 := 1/math.Pow(float64(k), s)
		t2 := RiemannZeta(s)
		p := t1/t2
		return p
	}
}

func Zeta_PMF_At(s float64, k int64) float64 {
	pmf := Zeta_PMF(s)
	return pmf(k)
}

func Zeta_CDF(s float64) func(k int64) float64 {
	return func(k int64) float64 {
		t1 := H(k, s)
		t2 := RiemannZeta(s)
		p := t1/t2
		return p
	}
}

func Zeta_CDF_At(s float64, k int64) float64 {
	pmf := Zeta_CDF(s)
	return pmf(k)
}

// Zeta distributed random variate
// Devroye 1986: 550. Called "Zipf distribution" there.
// Devroye, L. 1986: Non-Uniform Random Variate Generation. Springer-Verlag, New York. ISBN 0-387-96305-7.
func NextZeta(s float64) (k int64) {
	var x float64
	b := math.Pow(2.0, s - 1.0)
        for   {
                u := rand.Float64()
                v := rand.Float64()
                x = math.Floor(math.Pow(u, -1/(s - 1)))
                t := math.Pow(1 + 1.0/x, s - 1)
		delta := v*x*(t - 1.0)/(b - 1.0) 
		if delta <= (t/b) {
			break
		}
        }
	k = int64(x)
        return
}


func Zeta(s float64) func() int64 {
	return func() int64 { return NextZeta(s) }
}

func ZetaMean(s float64) float64 {
        if s <= 2 {
		panic("not defined")
	}
		t1 := RiemannZeta(s-1)
		t2 := RiemannZeta(s)
		return t1/t2
}

func ZetaMode() float64 {
		return 1
}

func ZetaVar(s float64) float64 {
        if s <= 3 {
		panic("not defined")
	}
		t1 := RiemannZeta(s)
		t2 := RiemannZeta(s-2)
		t3 := RiemannZeta((s-1)*(s-1))
		t4 := RiemannZeta(s*s)
		return (t1*t2-t3)/t4
}

