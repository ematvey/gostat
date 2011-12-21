// test for Quantile functions
// uses CDF and back using Qtl

package stat

import (
	"testing"
	"fmt"
	"math"
	"math/rand"
)

// test of Beta distribution
func Test_Quantiles_Beta(t *testing.T) {
	const nRuns int64 = 10
	var i int64
	fmt.Println("test for Quantile functions:")
	fmt.Println("\tBeta distribution")
	for i=0; i<nRuns; i++ {
		α := 6*rand.Float64()+0.05
		β := 6*rand.Float64()+0.05
		x := 0.9*rand.Float64()+0.05
		p := Beta_CDF_At(α, β, x)
		if p > 1e-15 {

		y := Beta_Qtl_For(α, β, p)
		if !check(x, y){
			t.Error()
			fmt.Println(α, β, p, x, y)
		}
		}
	}
}

// test of F distribution
func Test_Quantiles_F(t *testing.T) {
	const nRuns int64 = 10
	var i int64
	fmt.Println("\tF-distribution")
	for i=0; i<nRuns; i++ {
		df1 := 100*rand.Float64()+1.05
		df2 := 100*rand.Float64()+1.05
		x := 5*rand.Float64()
		p := F_CDF_At(math.Ceil(df1), math.Ceil(df2), x)
		if p > 1e-15 {
			y := F_Qtl_For(math.Ceil(df1), math.Ceil(df2), p)
			if !check(x, y){
				t.Error()
				fmt.Println(math.Ceil(df1), math.Ceil(df2), p, x, y)
			}
		}
	}
}

// test of Gamma distribution
func Test_Quantiles_Gamma(t *testing.T) {
	const nRuns int64 = 100
	var i int64
	fmt.Println("\tGamma distribution")
	for i=0; i<nRuns; i++ {
		k := 10*rand.Float64()+1.05
		θ := 3*rand.Float64()+0.05
		x := 25*rand.Float64()
		p := Gamma_CDF_At(math.Ceil(k), θ, x)
		if p > 1e-15 && p < 1- 1e-15{
			y := Gamma_Qtl_For(math.Ceil(k), θ, p)
			if !check(x, y){
				t.Error()
				fmt.Println(math.Ceil(k), θ, p, x, y)
			}
		}
	}
}

// test of Standard Normal distribution
func Test_Quantiles_Z(t *testing.T) {
	const nRuns int64 = 100
	var i int64
	fmt.Println("\tStandard Normal distribution")
	for i=0; i<nRuns; i++ {
		x := 6*rand.Float64()
		p := Z_CDF_At(x)
		if p > 1e-15 && p < 1- 1e-15{
			y := Z_Qtl_For(p)
			if !check(x, y){
				t.Error()
				fmt.Println(p, x, y)
			}
		}
	}
}



