// Empirical quantile from a sample from a PDF | PMF.
package bayes
import "sort"

// Empirical quantile from a sample from a PDF | PMF.
// returns α-quantile
func eQtl(x []float64, α float64) float64 {
	var (
		i int
		v float64
	)

	n := len(x)
	sort.Float64s(x)
	p := 0.0
	for i = 0; i < n && p < α*float64(n); i++ {
		p++
	}
	if i == n {
		v = x[i]
	} else {	// linear interpolation
		dp := α*float64(n) - p
		dx := (x[i+1]-x[i]) * dp
		v = x[i] + dx
	}
	return v
}

