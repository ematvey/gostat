// test of Next functions
// uses mean and variance of a big sample against computed values

package stat

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
)

// test of Normal distribution
func Test_Next_Normal(t *testing.T) {
	const nRuns int = 1
	const iter int = 100000	// need to be at least 1e5
	var n int
	var x, μ, σ, μ2, σ2, m, m2, delta float64
	fmt.Println("tests of Next functions:")
	fmt.Println("\tNormal distribution")

	for i := 0; i < nRuns; i++ {
		μ = 6 * rand.Float64()
		σ = 6*rand.Float64() + 1e-15
		μ2 = 0.0	// sample mean
		σ2 = 0.0	// sample variance unbiased
		m = 0.0
		m2 = 0.0

		for j := 0; j < iter; j++ {
        		n += 1
			x =  NextNormal(μ, σ)
			μ2 += x
			delta = x - m
        		m += delta/float64(n)
        		m2 += delta*(x - m)
		}
		σ2 = math.Sqrt(m2/float64(n - 1))
		μ2 /= float64(iter)

		if !check(μ, μ2) {
				t.Error()
				fmt.Println(μ, μ2)
		}
		if !check(σ, σ2) {
				t.Error()
				fmt.Println(σ, σ2)
		}
	}
}


// test of Gamma distribution
func Test_Next_Gamma(t *testing.T) {
	const nRuns int = 3
	const iter int = 10000	// need to be at least 1e5
	var n int
	var x, k, θ, μ, var1, μ2, var2, m, m2, delta float64
	fmt.Println("\tGamma distribution")

	for i := 0; i < nRuns; i++ {
		k = math.Ceil(6 * rand.Float64())
		θ = 6*rand.Float64() + 1e-15
		μ = GammaMean(k, θ)
		var1 = GammaVar(k, θ)
		μ2 = 0.0	// sample mean
		var2 = 0.0	// sample variance unbiased
		m = 0.0
		m2 = 0.0

		for j := 0; j < iter; j++ {
        		n += 1
			x =  NextGamma(k, θ)
			μ2 += x
			delta = x - m
        		m += delta/float64(n)
        		m2 += delta*(x - m)
		}
		var2 = math.Sqrt(m2/float64(n - 1))
		μ2 /= float64(iter)

		if !check(μ, μ2) {
				t.Error()
				fmt.Println(k, θ, μ, μ2)
		}
		if !check(var1, var2) {
				t.Error()
				fmt.Println(k, θ, var1, var2)
		}
	}
}



