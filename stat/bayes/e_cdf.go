// Empirical CDF
package bayes
import "sort"

// Empirical CDF, x is an array of observed values, y is a value where CDF is evaluated
// returns CDF value at y
func eCDF(x []float64, y float64) float64 {
	var (
		i int
		v float64
	)

	n := len(x)
	sort.Float64s(x)
	for i = 0; i < n && x[i] < y; i++ {
	}
	//return (2*float64(i)+1)/(2*float64(n))	// linear interpolation


	if i == n {
		v = 1
	} else {	// linear interpolation
		dx := y - x[i]
		dp := dx/(x[i+1]-x[i])
		v = (float64(i)+ dp)/float64(n)
	}
	return v

}


