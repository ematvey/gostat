// Empirical CDF
package bayes

// Empirical CDF, x is an array of observed values, y is a value where CDF is evaluated
func ecdf(x []float64, y float64) float64 {
		p := 0.0
		for i := 0; i < len(x) && x[i] < y; i++ {
			p += x[i]
		}
	return p
}


