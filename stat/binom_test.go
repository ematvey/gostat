package stat

import (
	"testing"
	"fmt"
)

// test against known values
func TestBinomial_PMF_CDF(t *testing.T) {
	var (
	α, prob float64
	i, n int64
	)


	// edit the following values:  >>>
	α=0.5
	n=20

k:=[]int64{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}
pmf:=[]float64{9.5367431640625E-007,1.9073486328125E-005,0.0001811981,0.0010871887,0.0046205521,0.0147857666,0.0369644165,0.073928833,0.1201343536,0.1601791382,0.176197052,0.1601791382,0.1201343536,0.073928833,0.0369644165,0.0147857666,0.0046205521,0.0010871887,0.0001811981,1.9073486328125E-005,9.5367431640625E-007}
cdf:=[]float64{9.5367431640625E-007,2.00271606445312E-005,0.0002012253,0.001288414,0.0059089661,0.0206947327,0.0576591492,0.1315879822,0.2517223358,0.411901474,0.588098526,0.7482776642,0.8684120178,0.9423408508,0.9793052673,0.9940910339,0.998711586,0.9997987747,0.9999799728,0.9999990463,1}

	// <<<

	fmt.Println("test of Binomial PMF")
	for i = 0; i < int64(len(k)); i++ {
		prob=Binomial_PMF_At(α, n, k[i])
			if !check(prob, pmf[i]){
				t.Error()
				fmt.Println(k[i], prob, pmf[i])

			}
	}
	fmt.Println("test of Binomial CDF")
	for i = 0; i < int64(len(k)); i++ {
		prob=Binomial_CDF_At(α, n, k[i])
			if !check(prob, cdf[i]){
				t.Error()
				fmt.Println(k[i], prob, cdf[i])
			}
	}
}





