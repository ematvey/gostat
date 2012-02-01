// Sample mean and variance (unbiased estimator)

package stat

import (
	"math"
)

// Sample mean and unbiased (Bessel correction) variance estimates
func SampleMeanVar(x []float64) (μ float64, σ float64) {
	var n int
	var m, m2, delta float64
		μ = 0.0	// sample mean
		σ = 0.0	// sample variance unbiased
		m = 0.0
		m2 = 0.0

		for j := 0; j < len(x); j++ {
        		n += 1
			μ += x[j]
			delta = x[j] - m
        		m += delta/float64(n)
        		m2 += delta*(x[j] - m)
		}
		σ = math.Sqrt(m2/float64(n - 1))
		μ /= float64(len(x))
	return
}

