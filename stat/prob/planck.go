// Planck distribution
// Johnson and Kotz, 1970
// x > 0.0 && a > 0.0 && b > 0.0

package prob

import (
	"math"
	. "go-fn.googlecode.com/hg/fn"
)


// ζ() waiting for better implementation
// Probability Density Function for the Planck distribution
func Planck_PDF(a, b float64) func(x float64) float64 {
	ζ := RiemannZeta
	return func(x float64) float64 {
	t1 := math.Pow(b, a+1)
	t2 := math.Pow(x, a)
	t3 := Γ(a+1) * ζ(a+1)
	t4 := math.Exp(b*x) - 1
		p := (t1*t2)/(t3*t4)
		return p
	}
}

// Planck_ distributed random variate
// Devroye 1986: 552.
// Devroye, L. 1986: Non-Uniform Random Variate Generation. Springer-Verlag, New York. ISBN 0-387-96305-7.

func NextPlanck(a, b float64) (x float64) {
	g := NextGamma(a+1, 1) // <<<  ??? OK ??
	z := float64(NextZeta(a+1))
        return g/(b*z)
}

func Planck(a, b float64) func() float64 {
	return func() float64 { return NextPlanck(a, b) }
}



