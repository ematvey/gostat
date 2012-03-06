// Credible interval for a sample from a posterior density
package bayes

// Credible interval for a sample from a posterior density
func ECrI(𝛩 []float64, width float64) (lo, hi float64) {
	α := (1-width)
	lo = eQtl(𝛩 , α/2)
	hi = eQtl(𝛩 , 1-α/2) 
	return 
}

