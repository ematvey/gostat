// Zipf distribution

package prob

import (
	"math"
	"math/rand"
)


// Zipf distributed random number
// Devroye 1986: 550.
// Devroye, L. 1986: Non-Uniform Random Variate Generation. Springer-Verlag, New York. ISBN 0-387-96305-7.
func NextZipf(a float64) (x int64) {
	b := math.Pow(2.0, a - 1.0)
        for   {
                u := rand.Float64()
                v := rand.Float64()
                x := math.Floor(math.Pow(u, -1.0/a - 1.0))
                t := math.Pow(1.0 + 1.0/x, a - 1.0)
		delta := v*x*(t - 1.0)/(b - 1.0) 
		if delta < (t/b) {
			break
		}
        }
        return
}



