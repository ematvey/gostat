package main

import (
	"fmt"
	"math"
	"math/rand"
	"gostat.googlecode.com/hg/stat"
)
func main() {  // just for testing purposes
	const z float64 = 1e2
	const m float64 = 5
	const tests = 1e2
	var λ, p, err, α, β, a, b, x, inv_cdf float64
	var count, n, k int64


	// test for Beta_CDF
	fmt.Println("test for Beta_CDF")
	fmt.Println("")
	α=2.0
	β=5.0
	cdf:=stat.Beta_CDF(α, β)
	fmt.Println("This should equal to 1.00:", cdf(0.999999999999999999))
	β=2.0
	cdf=stat.Beta_CDF(α, β)
	fmt.Println("This should equal to  0.50:", cdf(0.5))

	// test for Binomial
	fmt.Println("test for Binomial")
	fmt.Println("")
	p = 0.3897
	n = 25
	k = 9
	pmf := stat.Binomial_PMF(p, n)
	fmt.Println("p =", p, "  n=", n, "  k=", k, "  pmf =", pmf(k), " should be ")
	cdf2 := stat.Binomial_CDF(p, n)
	fmt.Println("p =", p, "  n=", n, "  k=", k, "  cdf =", cdf2(k), " should be ")


	// test for Poisson
	fmt.Println("test for Poisson_CDF")
	fmt.Println("# of tests: ", tests)
	fmt.Println("")
	for count = 0; count < tests;  {
		λ=z*rand.Float64()
		cdf3:=stat.Poisson_CDF(λ)
		k=(int64)(m*λ)
		p=cdf3(k)
		err=math.Abs(p-1)
			if (err > 1e-12)  {
				fmt.Println("λ =",λ , "  p =", p, "  k=", k)
			}
		count++
	}


	// test for InverseBeta CDF
	fmt.Println("test for InverseBeta CDF")
	fmt.Println("# of tests (except of those that did not converge: ", tests)
	fmt.Println("")
	for count := 0; count < tests;  {
		a=6*rand.Float64()+0.3
		b=6*rand.Float64()+0.3
		x=rand.Float64()
		p=stat.Beta_CDF_At(a, b, x)
		inv_cdf=stat.BetaInv_CDF_For(a, b, p)
		err=math.Abs(inv_cdf - x)
		if math.Abs(inv_cdf) < 2.0 && p < 1.00 {
			count++
			if (err > 1e-4)  {
				fmt.Println("α =",a , "  β =", b, "  p =", p, "  x =", x, "  err=", err, "  inv_cdf=",  inv_cdf)
			}
		}
	}
}

