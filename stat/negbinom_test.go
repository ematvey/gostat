package stat

import (
	"testing"
	"fmt"
)

// test against known values
func TestNegativeBinomial_PMF_CDF(t *testing.T) {
	var (
	ρ, prob float64
	i, n int64
	)


	// edit the following values:  >>>
	ρ=0.5
	n=20

k:=[]int64{10, 11, 12, 16, 25, 40}
pmf:=[]float64{0.0186544004827737808228, 0.025437818840146064758, 0.0328571826685220003128, 0.05907974191359244287, 0.04004139896255765052, 0.00121194851197753156874}
cdf:=[]float64{0.0493685733526945114136, 0.074806392192840576172, 0.1076635748613625764847, 0.30885965851484797895, 0.81435098276449480181, 0.9968911986703366647292}

	// <<<

	fmt.Println("test for NegativeBinomial PMF")
	for i = 0; i < int64(len(k)); i++ {
		prob=NegativeBinomial_PMF_At(ρ, n, k[i])
			if !check(prob, pmf[i]){
				t.Error()
				fmt.Println(k[i], prob, pmf[i])

			}
	}

	fmt.Println("test for NegativeBinomial CDF")
	for i = 0; i < int64(len(k)); i++ {
		prob=NegativeBinomial_CDF_At(ρ, n, k[i])
			if !check(prob, cdf[i]){
				t.Error()
				fmt.Println(k[i], prob, cdf[i])
			}
	}
}

