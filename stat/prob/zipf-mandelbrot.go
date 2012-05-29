// Zipf-Mandelbrot distribution
// For q == 0 it reduces to Zipf distribution.

package prob

import (
	"math"
	"math/rand"
)

// Check of Zipf-Mandelbrot function parameters
func ZipfMandelbrotCheckParams(n int64, q, s float64) bool {
	v:= false
	if n >0 && q >= 0 && s > 0 {
		v = true
	}
	return v
}

// Generalized harmonic number
func h(n int64, q, s float64) float64 {
	var i int64
	h := 0.0
	for i = 1; i <= n; i++ {
		h += math.Pow((float64(i) + q), -s)
	}
	return h
}

// Probability Mass Function for the Zipf-Mandelbrot distribution
func ZipfMandelbrotPMF(n int64, q, s float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := 1 / (math.Pow((float64(k)+q), s) * h(n, q, s))
		return p
	}
}

// Probability Mass Function for the Zipf-Mandelbrot distribution at k
func ZipfMandelbrotPMF_At(n int64, q, s float64, k int64) float64 {
	pmf := ZipfMandelbrotPMF(n, q, s)
	return pmf(k)
}

// Cumulative Distribution Function for the Zipf-Mandelbrot distribution
func ZipfMandelbrotCDF(n int64, q, s float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := h(k, q, s) / h(n, q, s)
		return p
	}
}

// Cumulative Distribution Function for the Zipf-Mandelbrot distribution at k
func ZipfMandelbrotCDF_At(n int64, q, s float64, k int64) float64 {
	cdf := ZipfMandelbrotCDF(n, q, s)
	return cdf(k)
}

// Quantile Function for the Zipf-Mandelbrot distribution
func ZipfMandelbrotQtl(n int64, q, s float64) func(p float64) int64 {
	return func(p float64) int64 {
		var k int64
		const kMax = 1e16
		cdf := ZipfMandelbrotCDF(n, q, s)
		if cdf(1) >= p {
			k = 1
		} else {
			for k=1; cdf(k) < p ; k++ {
				if k > kMax {
					panic("not found")
				}
			}
		}
		return k
	}
}

// Zipf-Mandelbrot distributed random variate
func NextZipfMandelbrot(n int64, q, s float64) (k int64) {
	qtl := ZipfMandelbrotQtl(n, q, s)
	p:= rand.Float64()
        return qtl(p)
}

// Zipf-Mandelbrot distribution function
func ZipfMandelbrot(n int64, q, s float64) func() int64 {
	return func() int64 { return NextZipfMandelbrot(n, q, s) }
}

// Mean of the Zipf-Mandelbrot distribution
func ZipfMandelbrotMean(n int64, q, s float64) float64 {
	return h(n, q, s-1)/h(n, q, s) - q
}


