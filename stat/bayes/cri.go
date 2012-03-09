// Bayesian credible interval. 

package bayes
/*
import (
	. "gostat.googlecode.com/hg/stat/prob"
	"math"
)
*/

// Bayesian credible interval for (analytical) quantile function 
func CrI(α float64, qtl func(𝛩 float64) float64) (hi, lo float64) {
	p := (1-α)
	lo = qtl(p/2)
	hi = qtl(1-p/2) 
	return
}

// Credible interval for a sample from a posterior density
func ECrI(𝛩 []float64, α float64) (lo, hi float64) {
	p := (1-α)
	lo = eQtl(𝛩 , p/2)
	hi = eQtl(𝛩 , 1-p/2) 
	return 
}

