package stat

import (
	"testing"
	"fmt"
)

// test against known values
func TestPoisson_PMF_CDF(t *testing.T) {
	var (
	α, prob float64
	i int64
	)


	// edit the following values:  >>>
	α=6
	k:=[]int64{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}
	pmf:=[]float64{0.0024787522,0.0148725131,0.0446175392,0.0892350784,0.1338526175,0.160623141,0.160623141,0.137676978,0.1032577335,0.068838489,0.0413030934,0.02252896,0.01126448,0.0051989908,0.0022281389,0.0008912556,0.0003342208,0.0001179603,3.93200983298978E-005,1.24168731568098E-005,3.72506194704295E-006}
	cdf:=[]float64{0.0024787522,0.0173512652,0.0619688044,0.1512038828,0.2850565003,0.4456796414,0.6063027824,0.7439797605,0.847237494,0.916075983,0.9573790764,0.9799080365,0.9911725165,0.9963715073,0.9985996462,0.9994909017,0.9998251226,0.9999430829,0.999982403,0.9999948198,0.9999985449}
	// <<<

	fmt.Println("test of Poisson PMF")
	for i = 0; i < int64(len(k)); i++ {
		prob=Poisson_PMF_At(α, k[i])
			if !check(prob, pmf[i]){
				t.Error()
				fmt.Println(k[i], prob, pmf[i])

			}
	}
	fmt.Println("test of Poisson CDF")
	for i = 0; i < int64(len(k)); i++ {
		prob=Poisson_CDF_At(α, k[i])
		if !check(prob, cdf[i]){
				t.Error()
				fmt.Println(k[i], prob, cdf[i])
		}

		cf:=Poisson_CDF_a(α)
		if !check(cf(k[i]), cdf[i]){
				t.Error()
				fmt.Println(k[i], prob, cdf[i])
		}
		
	}
}

