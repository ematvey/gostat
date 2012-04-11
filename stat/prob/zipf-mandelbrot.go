// Zipf-Mandelbrot distribution

package prob

import (
	"math"
)

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

/*
func NextZipf-Mandelbrot(n int64, q, s float64) (k int64) {
        return 
}

func Zipf-Mandelbrot(n int64, q, s float64) func() int64 {
	return func() int64 { return NextZipfMandelbrot(n, q, s) }
}
*/

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
