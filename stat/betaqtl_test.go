// test for Beta_Qtl_For(α, β, p)
package stat

import (
	"testing"
	"fmt"
	"math"
	"math/rand"
)

// test against known values
func TestBeta_Qtl_For(t *testing.T) {
	fmt.Println("test for Beta_Qtl_For")
	var x, y, err, α, β, p float64
	var count, tests int64

	α=10.001
	β=5.0001
	p=0.01
	x=Beta_Qtl_For(α, β, p)
	y=0.3726
	if !check(x, y){
		t.Error()
	}

	p=0.5
	x=Beta_Qtl_For(α, β, p)
	y=0.6742
	if !check(x, y){
		t.Error()
	}

	p=0.99
	x=Beta_Qtl_For(α, β, p)
	y=0.8981
	if !check(x, y){
		t.Error()
	}

	for count = 0; count < tests;  {
		α=6*rand.Float64()+0.3
		β=6*rand.Float64()+0.3
		x=rand.Float64()
		p=Beta_CDF_At(α, β, x)
		inv_cdf:=Beta_Qtl_For(α, β, p)
		err=math.Abs(inv_cdf - x)
		if math.Abs(inv_cdf) < 2.0 && p < 1.00 {
			count++
			if !check(inv_cdf, x){
				t.Error()
				fmt.Println("α =",α , "  β =", β, "  p =", p, "  x =", x, "  err=", err, "  inv_cdf=",  inv_cdf)
			}
		}
	}
}


