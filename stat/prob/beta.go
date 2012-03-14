// Beta distribution

package prob

import (
	"fmt"
	. "go-fn.googlecode.com/hg/fn"
	"math"
)

func bisect(x, p, a, b, xtol, ptol float64) float64 {

	var x0, x1, px float64

	cdf := Beta_PDF(a, b)

	for math.Abs(x1-x0) > xtol {
		px = cdf(x)
		switch {
		case math.Abs(px-p) < ptol:
			return x
		case px < p:
			x0 = x
		case px > p:
			x1 = x
		}
		x = 0.5 * (x0 + x1)
	}
	return x
}

func betaContinuedFraction(α, β, x float64) float64 {

	var aa, del, res, qab, qap, qam, c, d, m2, m, acc float64
	var i int64
	const eps = 2.2204460492503131e-16
	const maxIter = 1000000000

	acc = 1e-16
	qab = α + β
	qap = α + 1.0
	qam = α - 1.0
	c = 1.0
	d = 1.0 - qab*x/qap

	if math.Abs(d) < eps {
		d = eps
	}
	d = 1.0 / d
	res = d

	for i = 1; i <= maxIter; i++ {
		m = (float64)(i)
		m2 = 2 * m
		aa = m * (β - m) * x / ((qam + m2) * (α + m2))
		d = 1.0 + aa*d
		if math.Abs(d) < eps {
			d = eps
		}
		c = 1.0 + aa/c
		if math.Abs(c) < eps {
			c = eps
		}
		d = 1.0 / d
		res *= d * c
		aa = -(α + m) * (qab + m) * x / ((α + m2) * (qap + m2))
		d = 1.0 + aa*d
		if math.Abs(d) < eps {
			d = eps
		}
		c = 1.0 + aa/c
		if math.Abs(c) < eps {
			c = eps
		}
		d = 1.0 / d
		del = d * c
		res *= del
		if math.Abs(del-1.0) < acc {
			return res
		}
	}

	panic(fmt.Sprintf("betaContinuedFraction(): α or β too big, or maxIter too small"))
	return -1.00
}

func Beta_PDF(α float64, β float64) func(x float64) float64 {
	dα := []float64{α, β}
	dirPDF := Dirichlet_PDF(dα)
	return func(x float64) float64 {
		if 0 > x || x > 1 {
			return 0
		}
		dx := []float64{x, 1 - x}
		return dirPDF(dx)
	}
}
func Beta_LnPDF(α float64, β float64) func(x float64) float64 {
	dα := []float64{α, β}
	dirLnPDF := Dirichlet_LnPDF(dα)
	return func(x float64) float64 {
		if 0 > x || x > 1 {
			return negInf
		}
		dx := []float64{x, 1 - x}
		return dirLnPDF(dx)
	}
}
func NextBeta(α float64, β float64) float64 {
	dα := []float64{α, β}
	return NextDirichlet(dα)[0]
}
func Beta(α float64, β float64) func() float64 {
	return func() float64 { return NextBeta(α, β) }
}

// Value of PDF of Beta distribution(α, β) at x
func Beta_PDF_At(α, β, x float64) float64 {
	pdf := Beta_PDF(α, β)
	return pdf(x)
}

// CDF of Beta-distribution
func Beta_CDF(α float64, β float64) func(x float64) float64 {
	return func(x float64) float64 {
		//func Beta_CDF(α , β , x float64) float64 {
		var y, res float64
		y = math.Exp(LnΓ(α+β) - LnΓ(α) - LnΓ(β) + α*math.Log(x) + β*math.Log(1.0-x))
		switch {
		case x == 0:
			res = 0.0
		case x == 1.0:
			res = 1.0
		case x < (α+1.0)/(α+β+2.0):
			res = y * betaContinuedFraction(α, β, x) / α
		default:
			res = 1.0 - y*betaContinuedFraction(β, α, 1.0-x)/β

		}
		return res
	}
}

// Value of CDF of Beta distribution(α, β) at x
func Beta_CDF_At(α, β, x float64) float64 {
	cdf := Beta_CDF(α, β)
	res := cdf(x)
	return res
}

// Inverse of the cumulative beta probability density function (quantile).
// p: Probability associated with the beta distribution
// α: Parameter of the distribution
// β: Parameter of the distribution
func Beta_Qtl(α, β float64) func(p float64) float64 {
	return func(p float64) float64 {
		var x float64 = 0
		var a float64 = 0
		var b float64 = 1
		var precision float64 = 1e-9
		if p < 0.0 {
			panic(fmt.Sprintf("p < 0"))
		}
		if p > 1.0 {
			panic(fmt.Sprintf("p > 1.0"))
		}
		if α < 0.0 {
			panic(fmt.Sprintf("α < 0.0"))
		}
		if β < 0.0 {
			panic(fmt.Sprintf("β < 0.0"))
		}

		for (b - a) > precision {
			x = (a + b) / 2
			if BetaIncReg(α, β, x) > p {
				b = x
			} else {
				a = x
			}
		}

		return x
	}
}

// Inverse of the cumulative beta probability density function (quantile) for a given probability.
// p: Probability associated with the beta distribution
// α: Parameter of the distribution
// β: Parameter of the distribution
func Beta_Qtl_For(α, β, p float64) float64 {
	cdf := Beta_Qtl(α, β)
	return cdf(p)
}
