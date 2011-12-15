package stat

import (
	"time"
	"math"
	"math/rand"
	"testing"
	"fmt"
	. "go-fn.googlecode.com/hg/fn"
)

var Seed func(int64) = rand.Seed

func XTestDir(t *testing.T) {
	α := []float64{4, 5, 6}
	dgen := Dirichlet(α)
	counts := [3]int{0, 0, 0}
	const total = 150000
	for i := 0; i < total; i++ {
		θ := dgen()
		v := NextChoice(θ)
		counts[v]++
	}
	fmt.Printf("%v\n", counts)
}

func TestNullWeights(t *testing.T) {
	n := int64(10)
	weights := make([]float64, n)
	m := NextChoice(weights)
	if n != m {
		t.Error()
	}
}

func TestLnGamma(t *testing.T) {
	acc := 0.0000001
	check := func(x, y float64) bool {
		if false {
			return x == y
		}
		return math.Abs(x-y) < acc
	}
	for i := 0; i < 100; i++ {
		x := NextGamma(10, 10)
		g1 := LnΓ(x)
		g2, _ := math.Lgamma(x)
		if !check(g1, g2) {
			t.Error(fmt.Sprintf("For %v: %v vs %v", x, g1, g2))
		}
	}
	//var start int64
	Seed(10)
	start := time.Now()
	for i := 0; i < 1e6; i++ {
		x := NextGamma(10, 10)
		math.Lgamma(x)
	}
	now:=time.Now()
	duration2 := float64(now.Sub(start)) / 1e9

	//duration2 := float64(time.Now()-start) / 1e9
	Seed(10)
	start = time.Now()
	for i := 0; i < 1e6; i++ {
		x := NextGamma(10, 10)
		LnΓ(x)
	}
	now=time.Now()
	duration1 := float64(now.Sub(start)) / 1e9
	fmt.Printf("Mine was %f\nTheirs was %f\n", duration1, duration2)
}

func XTestGen(t *testing.T) {
	fmt.Printf("NextUniform => %f\n", NextUniform())
	fmt.Printf("NextExp => %f\n", NextExp(1.5))
	fmt.Printf("NextGamma => %f\n", NextGamma(.3, 1))
	fmt.Printf("NextNormal => %f\n", NextNormal(0, 1))
	fmt.Printf("NextRange => %d\n", NextRange(10))
	fmt.Printf("NextChoice => %d\n", NextChoice([]float64{.3, .3, .4}))
	fmt.Printf("NextMultinomial => %v\n",
		NextMultinomial([]float64{.3, .3, .4}, 100))
	fmt.Printf("NextDirichlet => %v\n", NextDirichlet([]float64{.3, .3, .4}))
	fmt.Printf("NextBernoulli => %d\n", NextBernoulli(.5))
	fmt.Printf("NextGeometric => %d\n", NextGeometric(.5))
	fmt.Printf("NextBinomial => %d\n", NextBinomial(.5, 10))
	fmt.Printf("NextPoisson => %d\n", NextPoisson(1.5))
	fmt.Printf("NextXsquare => %f\n", NextXsquare(3))
	fmt.Printf("NextNegativeBinomial => %d\n", NextNegativeBinomial(.5, 10))
	fmt.Printf("NextStudentsT => %f\n", NextStudentsT(7))
	fmt.Printf("NextF => %f\n", NextF(7, 3))
/*	fmt.Printf("NextWishart => %v\n",

		NextWishart(100, matrix.MakeDenseMatrixStacked([][]float64{[]float64{1, 0}, []float64{0, 1}})))
	fmt.Printf("NextInverseWishart => %v\n",
		NextInverseWishart(100, matrix.MakeDenseMatrixStacked([][]float64{[]float64{1, 0}, []float64{0, 1}})))
*/
}

// test for Beta_CDF
func TestBeta_CDF(t *testing.T) {
	var x, y, acc, α, β float64
	acc = 1e-5
	check := func(x, y, acc float64) bool {
		if math.Abs(x-y) > acc  {
			return false
		}
		return true
	}

	α=2.0
	β=5.0
	cdf:=Beta_CDF(α, β)
	x=cdf(0.999999999999999999)
	y=1.00

	if !check(x, y, acc){
		t.Error()
	}

	β=2.0
	cdf=Beta_CDF(α, β)
	x=cdf(0.5)
	y=0.50

	if !check(x, y, acc){
		t.Error()
	}

	α=0.5
	β=2.0
	cdf=Beta_CDF(α, β)
	x=cdf(0.1)
	y=0.4585302607

	if !check(x, y, acc){
		t.Error()
	}

	x=cdf(0.6)
	y=0.9295160031

	if !check(x, y, acc){
		t.Error()
	}

	x=cdf(0.7)
	y=0.9621590305

	if !check(x, y, acc){
		t.Error()
	}

}

// test for Binomial
func TestBinomial(t *testing.T) {
	var x, y, acc, p float64
	var n, k int64
	acc = 1e-5
	check := func(x, y, acc float64) bool {
		if math.Abs(x-y) > acc  {
			return false
		}
		return true
	}
	p = 0.3897
	n = 25
	k = 9

	pmf := Binomial_PMF(p, n)
	x = pmf(k)
	y = 0.1568666506

	if !check(x, y, acc){
		t.Error()
	}


	cdf := Binomial_CDF(p, n)
	x = cdf(k)
	y = 0.4665806719

	if !check(x, y, acc){
		t.Error()
	}

}


// test for Incomplete Beta
func TestIncompleteBeta(t *testing.T) {
	var x, y, acc float64
	acc = 1e-5
	check := func(x, y, acc float64) bool {
		if math.Abs(x-y) > acc  {
			return false
		}
		return true
	}

	x = IB(0.5, 2, 0.5)
	y = 1.17851130

	if !check(x, y, acc){
		t.Error()
	}

	x = IB(2.0, 6, 0.75)
	y = 0.02377755


	if !check(x, y, acc){
		t.Error()
	}

	x = IB(1.456, 3.25895, 0.99)
	y = 0.14450022


	if !check(x, y, acc){
		t.Error()
	}
}

// test for Incomplete Gamma
func TestIncompleteGamma(t *testing.T) {
	var x, y, z, acc float64
	acc = 1e-05
	check := func(x, y, acc float64) bool {
		if x/y > 1.00 {
			z = y/x
		} else {
			z = x/y
		}
		if 1-z > acc  {
			return false
		}
		return true
	}
	x = IΓ(3, 5.5)
	y = 0.17675286

	if !check(x, y, acc){
		t.Error()
	}

	x = IΓ(6, 3.96545)
	y = 94.86080842

	fmt.Println(x, y)

	if !check(x, y, acc){
		t.Error()
	}
/*
	x = IΓ(18.3542, 3.96545)
	y = 9.838284e+14
	fmt.Println(x, y)

	if !check(x, y, acc){
		t.Error()
	}
//FAILED, non-integer s
*/

	x = IΓ(18, 3.96545)
	y = 3.556874e+14

	if !check(x, y, acc){
		t.Error()
	}

	
}
// test for Poisson
func TestPoisson(t *testing.T) {
	var x, y, acc, λ float64
	var k int64
	acc = 1e-5
	check := func(x, y, acc float64) bool {
		var r float64
		if x/y > 1.00 {
			r = y/x
		} else {
			r = x/y
		}
		if 1-r > acc  {
			return false
		}
		return true
	}

	λ = 1.77
	k = 1

	pmf := Poisson_PMF(λ)
	x = pmf(k)
	y = 0.3014899390220975

	if !check(x, y, acc){
		t.Error()
		fmt.Println("k: ", k, "λ: ", λ, "prob: ", x, "err: ", 1-x/y)

	}

	fmt.Println("")
	fmt.Println("test for Poisson_CDF")
	fmt.Println("")


	λ = 0.15164076846159652
	k = 1
	cdf:=Poisson_CDF(λ)
	x=cdf(k)
	y=0.989601355904656
	if !check(x, y, acc){
		t.Error()
		fmt.Println("k: ", k, "λ: ", λ, "prob: ", x, "err: ", 1-x/y)
	}


}

// test for NegativeBinomial_CDF
func TestNegativeBinomial_CDF(t *testing.T) {
	fmt.Println("")
	fmt.Println("test for NegativeBinomial_CDF")
	fmt.Println("")
	var x, y, acc, ρ float64
	var r int64
	acc = 1e-5
	check := func(x, y, acc float64) bool {
		if math.Abs(x-y) > acc  {
			return false
		}
		return true
	}

	ρ = 0.5
	r = 20
	cdf5:=NegativeBinomial_CDF(ρ, r)
	x = cdf5(10)
	y = 0.0493685733526945114136
	if !check(x, y, acc){
		t.Error()
	}

	ρ = 0.5
	r = 10
	cdf5=NegativeBinomial_CDF(ρ, r)
	x = cdf5(20)
	y = 0.9786130273714661598206
	if !check(x, y, acc){
		t.Error()
	}
}

// test for BetaInv_CDF_For(α, β, p)
func TestBetaInv_CDF_For(t *testing.T) {
	fmt.Println("")
	fmt.Println("test for BetaInv_CDF_For(α, β, p)")
	fmt.Println("")
	var x, y, z, acc, err, α, β, p float64
	var count, tests int64

	acc = 1e-4

	check := func(x, y, acc float64) bool {
		if x/y > 1.00 {
			z = y/x
		} else {
			z = x/y
		}
		if 1-z > acc  {
			return false
		}
		return true
	}

	α=10.001
	β=5.0001
	p=0.01
	x=BetaInv_CDF_For(α, β, p)
	y=0.3726
	fmt.Println(x, " == ", y)

	if !check(x, y, acc){
		t.Error()
	}

	p=0.5
	x=BetaInv_CDF_For(α, β, p)
	y=0.6742
	fmt.Println(x, " == ", y)

	if !check(x, y, acc){
		t.Error()
	}

	p=0.99
	x=BetaInv_CDF_For(α, β, p)
	y=0.8981
	fmt.Println(x, " == ", y)

	if !check(x, y, acc){
		t.Error()
	}

	for count = 0; count < tests;  {
		α=6*rand.Float64()+0.3
		β=6*rand.Float64()+0.3
		x=rand.Float64()
		p=Beta_CDF_At(α, β, x)
		inv_cdf:=BetaInv_CDF_For(α, β, p)
		err=math.Abs(inv_cdf - x)
		if math.Abs(inv_cdf) < 2.0 && p < 1.00 {


			count++
			if !check(inv_cdf, x, acc){
				t.Error()
				fmt.Println("α =",α , "  β =", β, "  p =", p, "  x =", x, "  err=", err, "  inv_cdf=",  inv_cdf)
			}
		}
	}
}

// test for F_InvCDF_For
func TestF_InvCDF_For(t *testing.T) {
	fmt.Println("")
	fmt.Println("test for F_InvCDF_For")
	fmt.Println("")
	var acc, df1, df2, x, y, z, p float64
	acc = 1e-4

	check := func(x, y, acc float64) bool {
		if x/y > 1.00 {
			z = y/x
		} else {
			z = x/y
		}
		if 1-z > acc  {
			return false
		}
		return true
	}
	df1=3
	df2=3
	x=0.46
	cdf:=F_CDF(df1, df2)
	p=cdf(x)
	y = F_InvCDF_For(df1, df2, p)

	if !check(x, y, acc){
		t.Error()
	}

fmt.Println(F_InvCDF_For(3, 3, 0.1), " = 0.19")
fmt.Println(F_InvCDF_For(3, 3, 0.27002231), " = 0.46")
}


/*
// test for Binomial p confidence interval  // failed due to some unknown bug
func TestBinomP_CI(t *testing.T) {
	fmt.Println("")
	fmt.Println("test for Binomial p confidence interval")
	fmt.Println("")
	var n int64
	var k, nn, p, alpha, low, high, low2, high2 float64
	n=30
	nn=float64(n)
	p=0.1
	k=nn*p

	alpha=0.1
	low2=0.04
	high2=0.21
	low, high =  Binom_p_ConfI(n, p, alpha) 
	fmt.Println(low, " = ", low2, "\t", high, " = ",  high2)
}
*/

