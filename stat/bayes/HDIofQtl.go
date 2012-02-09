// Highest density iterval (HDI) limits, for quantile function
// Kruschke 2012: Chapter 23.3.3, p. 629 and further.

package bayes
import (
	"fmt"
	"math"
	. "gostat.googlecode.com/hg/stat"
)
// Auxiliary function
func aux(α, β, credMass float64) {
  func (x float64) { iw(α, β, credMass, x float64) }
}

  foo := aux(α, β, credMass)
  f:= foo( lowTailPr )


// Interval width
func iw(α, β, credMass, lowTailPr float64) {
      Beta_Qtl_For(α, β, credMass + lowTailPr) - Beta_Qtl_For(α, β, lowTailPr)
}

// HDI of Beta Distribution
func HDIofBetaQtl(α, β, credMass, tol float64) (lo, hi float64) {

fmin(ax,bx,f,tol)
	f := iw(α, β, credMass, lowTailPr)
	min = fmin(0 , 1.0 - credMass, f, tol)
	lo = Beta_Qtl_For(α, β, min)
	hi = Beta_Qtl_For(α, β, credMass + min)
	return	lo, hi
}

