package stat

import (
	"testing"
	"fmt"
)

// test against known values
func TestBeta_PDF_CDF(t *testing.T) {
	var (
		α, β, prob float64
		i int64
	)


	// edit the following values:  >>>
	α=0.5
	β=2.0

	x:=[]float64{0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1}
	pdf:=[]float64{3.1863968679,2.1345374206,1.6460179221,1.3416407865,1.125,0.9585144756,0.8240253984,0.7115124735,0.6149186938,0.5303300859,0.4550849072,0.3872983346,0.3255911783,0.2689264371,0.2165063509,0.1677050983,0.1220233825,0.0790569415,0.0384741882,0}
	cdf:=[]float64{0.3298200267,0.4585302607,0.5519001268,0.6260990337,0.6875,0.7394254526,0.7838805713,0.8221921916,0.8552960014,0.8838834765,0.9084843147,0.9295160031,0.9473152854,0.9621590305,0.9742785793,0.9838699101,0.9911010292,0.996117463,0.9990464203,1}

	// <<<

	fmt.Println("test of Beta PDF")
	for i = 0; i < int64(len(x)); i++ {
		prob=Beta_PDF_At(α, β, x[i])
			if !check(prob, pdf[i]){
				t.Error()
				fmt.Println(x[i], prob, pdf[i])

			}
	}

	fmt.Println("test of Beta CDF")
	for i = 0; i < int64(len(x)); i++ {
		prob=Beta_CDF_At(α, β, x[i])
			if !check(prob, cdf[i]){
				t.Error()
				fmt.Println(x[i], prob, cdf[i])
			}
	}
}

