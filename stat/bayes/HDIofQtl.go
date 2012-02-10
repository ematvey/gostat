// Highest density iterval (HDI) limits, for quantile function
// Kruschke 2012: Chapter 23.3.3, p. 629 and further.

package bayes
import (
	. "gostat.googlecode.com/hg/stat/pdf"
)
// Interval width
func iw(α, β, credMass, lowTailPr float64) float64{
      return Beta_Qtl_For(α, β, credMass + lowTailPr) - Beta_Qtl_For(α, β, lowTailPr)
}

// Interval width for fixed α, β, credMass
func iwFix(α, β, credMass float64) func (x float64) float64{
	return func (x float64) float64 { return iw(α, β, credMass, x)}
}

// HDI of Beta Distribution
func HDIofBetaQtl(α, β, credMass, tol float64) (lo, hi float64) {
	f := iwFix(α, β, credMass)
// func fmin(f func(float64) float64, ax, bx, tol float64) float64 {

	min := fmin(f, 0 , 1 - credMass, tol)
	lo = Beta_Qtl_For(α, β, min)
	hi = Beta_Qtl_For(α, β, credMass + min)
	return	lo, hi
}

