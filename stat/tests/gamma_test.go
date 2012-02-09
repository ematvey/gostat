package stat

import (
	"testing"
	"fmt"
)

// test against known values

func TestGamma_PDF_CDF(t *testing.T) {
	var (
		α, λ, prob float64
		i, k int64
	)


	// edit the following values:  >>>
	α=9.0
	k=9
	λ=0.5

	x:=[]float64{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}
	pdf:=[]float64{0,0.0017185433,0.0595403626,0.2065154671,0.2791730639,0.2251980643,0.1310465699,0.0608710805,0.0239749452,0.0083250881,0.0026173379,0.0007592981,0.0002061256,5.29224238287655E-005,1.29577092147652E-005,3.045404975127E-006,6.90694309390341E-007,1.51819803105039E-007,3.24584763738617E-008,6.76998852924674E-009,1.38105230394241E-009}
	cdf:=[]float64{0,0.0002374473,0.0213634345,0.152762506,0.4074526586,0.6671803212,0.8449722182,0.9379448041,0.9780127465,0.9929439909,0.997912741,0.9994230988,0.9998494373,0.9999625849,0.9999910876,0.9999979539,0.9999995453,0.9999999018,0.9999999793,0.9999999957,0.9999999991}
	// <<<

	fmt.Println("test of Gamma PDF")
	for i = 0; i < int64(len(x)); i++ {
		prob=Gamma_PDF_At(α, λ, x[i])
			if !check(prob, pdf[i]){
				t.Error()
				fmt.Println(x[i], prob, pdf[i])

			}
	}

/*
	fmt.Println("test of Gamma CDF")
	for i = 0; i < int64(len(x)); i++ {
		prob=Gamma_CDF_At(α, λ, x[i])
			if !check(prob, cdf[i]){
				t.Error()
				fmt.Println(x[i], prob, cdf[i])
			}
	}
*/
	fmt.Println("test of Gamma CDF with integer k")
	for i = 0; i < int64(len(x)); i++ {
		fn:=Gamma_CDFint(k, λ)
		prob=fn(x[i])
			if !check(prob, cdf[i]){
				t.Error()
				fmt.Println(x[i], prob, cdf[i])
			}
	}


}

