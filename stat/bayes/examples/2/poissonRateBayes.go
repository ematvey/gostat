// Summary of posterior of binomial parameter
package main

import (
	"fmt"
//	"math"
	"gostat.googlecode.com/hg/stat/bayes"
)


func main() {
	var x, n int64
	pr := []float64{0.005,0.01,0.025,0.05,0.5,0.95,0.975,0.99,0.995}
	r :=25.0
	v :=4.5
	x = 19 
	n = 5

	if r<0 || v<0 {
		panic("Shape parameter r and rate parameter v must be greater than or equal to zero")
	}
	fmt.Println("\nProb.\t\tQuantile \n")
	for i:=0; i< 9; i++ {
			qf := bayes.PoissonLambda_Qtl_GPri(x, n, r, v)
			qtl := qf(pr[i])
			fmt.Println(pr[i],"\t\t", qtl)
	}
	fmt.Println("\n")
}


