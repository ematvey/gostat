// Summary of posterior of binomial parameter
package main

import (
	"fmt"
	"math"
	"gostat.googlecode.com/hg/stat/bayes"
)


func main(){
	var k, n int64
	pr:= []float64{0.005,0.01,0.025,0.05,0.5,0.95,0.975,0.99,0.995}


	//fmt.Scanf("%v %v %v %v", k, n, a, b)
	k = 10
	n = 20
	a := 0.5
	b := 0.5

/*
  prior = dbeta(pi,a,b)
  likelihood = dbinom(k,n,prob=pi)
  posterior = dbeta(pi,a+k,b+n-k)
*/

	//* posterior summary
	m1 := (a+float64(k))/(a+b+float64(n))
	v1 := m1*(1-m1)/(a+b+float64(n)+1)
	s1 := math.Sqrt(v1)
	fmt.Println("Posterior Mean           : ", m1)
	fmt.Println("Posterior Variance       : ", v1)
	fmt.Println("Posterior Std. Deviation : ", s1)

	fmt.Println("Posterior Mean           : ", bayes.BinomPiPostMean(a,  b , n, k))
	fmt.Println("Posterior Variance       : ", bayes.BinomPiPostVar(a,  b , n, k))


	fmt.Println("\nProb.\t\tQuantile \n")
	for i:=0; i< 9; i++ {
			qf := bayes.BinomPi_Qtl_BPri(k, n, a, b)
			qtl := qf(pr[i])
			fmt.Println(pr[i],"\t\t", qtl)
	}
	fmt.Println("\n")
}


