// Generalized Pareto Distribution
// params:
// shape1 > 0.0
// shape2 > 0.0
// scale > 0.0
// support: x >= θ
// Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2008), Loss Models, From Data to Decisions, Third Edition, Wiley.

package prob

import (
	"math"
	"math/rand"
	. "go-fn.googlecode.com/hg/fn"
)

func ParetoG_PDF(shape1, shape2, scale float64) func(x float64) (p float64) {
	// We work with the density expressed as
	// u^shape2 * (1 - u)^shape1 / (x * B(shape1, shape2))
	// with u = v/(1 + v) = 1/(1 + 1/v), v = x/scale.
	return func(x float64) (p float64) {
		if x < 0 {
			p = 0
		} else if x == 0 {
			if shape2 < 1 {
				p = math.Inf(1)
			} else if shape2 > 1 {
				p = 0
			} else {
				p = 1 / (scale * B(shape2, shape1))
			}
		}
		tmp := math.Log(x) - math.Log(scale)
		logu := - math.Log1p(math.Exp(-tmp))
		log1mu := - math.Log1p(math.Exp(tmp))
		p =  math.Exp(shape2 * logu + shape1 * log1mu - math.Log(x) - LnB(shape2, shape1))
		return
	}
}

func ParetoG_PDF_At(shape1, shape2, scale, x float64) float64 {
	pdf := ParetoG_PDF(shape1, shape2, scale)
	return pdf(x)
}

func ParetoG_CDF(shape1, shape2, scale float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		u := math.Exp(-math.Log1p(math.Exp(math.Log(scale) - math.Log(x))))
		cdf := Beta_CDF(shape2, shape1)
		return cdf(u)
	}
}

func ParetoG_CDF_At(shape1, shape2, scale, x float64) float64 {
	cdf := ParetoG_CDF(shape1, shape2, scale)
	return cdf(x)
}

func ParetoG_Qtl(shape1, shape2, scale float64) func(p float64) float64 {
	return func(p float64) float64 {
		if p < 0 || p > 1 {
			panic("bad param")
		}
		qtl := Beta_Qtl(shape2, shape1)
		return scale / (1.0 / qtl(p) - 1.0)
	}
}

func NextParetoG(shape1, shape2, scale float64) float64 {
	qtl := ParetoG_Qtl(shape1, shape2, scale)
	p := rand.Float64()
	return qtl(p)
}

func ParetoG_Moment(shape1, shape2, scale float64, order int) (x float64) {
	o := float64(order)
	if o  <= -shape2 || o  >= shape1 {
		x = math.Inf(1)
	} else {
		x = math.Pow(scale, o) * B(shape1 - o, shape2 + o) / B(shape1, shape2)
	}
	return
}

func ParetoG_Mean(shape1, shape2, scale float64) float64 {
	return ParetoG_Moment(shape1, shape2, scale, 1)
}

func ParetoG_Var(shape1, shape2, scale float64) float64 {
	return ParetoG_Moment(shape1, shape2, scale, 2)
}

func ParetoG_Skew(shape1, shape2, scale float64) float64 {
	return ParetoG_Moment(shape1, shape2, scale, 3)
}

