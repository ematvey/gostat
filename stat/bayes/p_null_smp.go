// Probability of a one sided null hypothesis from a sample from a posterior density.
package bayes

// Lower tail probability of a one sided null hypothesis from a sample from a posterior density.
func PNullSmpLowT(𝛩 []float64, 𝛩0 float64) float64 {
	return eCDF(𝛩 , 𝛩0)
}

// Upper tail probability of a one sided null hypothesis from a sample from a posterior density.
func PNullSmpUppT(𝛩 []float64, 𝛩0 float64) float64 {
	return 1-eCDF(𝛩, 𝛩0)
}


