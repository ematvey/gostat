package prob

import (
	"testing"
	"fmt"
)

// test against known values
func TestHypergeometric_CDF(t *testing.T) {
	var (
	size, m, n, k int64
	)


	// edit the following values:  >>>
	size = 200
	m = 60
	n = 50
	k = 14
	known := 0.4339792084682835697
		
	prob :=	Hypergeometric_CDF_At(size, m, n, k)
	fmt.Println("test of Hypergeometric CDF")
	if !check(prob, known){
				t.Error()
				fmt.Println(prob, known)
	}


	size = 10
	m = 5
	n = 6
	k = 2
	known = 0.238095238095238095238
		
	prob =	Hypergeometric_PMF_At(size, m, n, k)
	fmt.Println("test of Hypergeometric PMF")
	if !check(prob, known){
				t.Error()
				fmt.Println(prob, known)
	}

	known = 0.26190476190476190476
			
	prob =	Hypergeometric_CDF_At(size, m, n, k)
	fmt.Println("2nd test of Hypergeometric CDF")
	if !check(prob, known){
				t.Error()
				fmt.Println(prob, known)
	}

	size = 100
	m = 49
	n = 55
	k = 27
	known = 0.15916742856360957427
		
		
	prob =	Hypergeometric_PMF_At(size, m, n, k)
	fmt.Println("test of Hypergeometric PMF")
	if !check(prob, known){
				t.Error()
				fmt.Println(prob, known)
	}

	known = 0.5873860392113934918
		
			
	prob =	Hypergeometric_CDF_At(size, m, n, k)
	fmt.Println("2nd test of Hypergeometric CDF")
	if !check(prob, known){
				t.Error()
				fmt.Println(prob, known)
	}
		
}

