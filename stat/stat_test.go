package stat

import (
	"testing";
	"fmt";
)

func TestDir(t *testing.T) {
	α := []float64{4,5,6}
	dgen := Dirichlet(α)
	counts := [3]int{0,0,0}
	const total = 150000
	for i:=0;i<total;i++ {
		θ := dgen()
		v := NextChoice(θ)
		counts[v]++
	}
	fmt.Printf("%v\n", counts)
}

func TestGen(t *testing.T) {
	fmt.Printf("NextUniform => %f\n", NextUniform());
	fmt.Printf("NextExp => %f\n", NextExp(1.5));
	fmt.Printf("NextGamma => %f\n", NextGamma(.3, 1));
	fmt.Printf("NextNormal => %f\n", NextNormal(0, 1));
	fmt.Printf("NextRange => %d\n", NextRange(10));
	fmt.Printf("NextChoice => %d\n", NextChoice([]float64{.3, .3, .4}));
	fmt.Printf("NextMultinomial => %v\n",
		NextMultinomial([]float64{.3, .3, .4}, 100));
	fmt.Printf("NextDirichlet => %v\n", NextDirichlet([]float64{.3, .3, .4}));
	fmt.Printf("NextBernoulli => %d\n", NextBernoulli(.5));
	fmt.Printf("NextGeometric => %d\n", NextGeometric(.5));
	fmt.Printf("NextBinomial => %d\n", NextBinomial(.5, 10));
	fmt.Printf("NextPoisson => %d\n", NextPoisson(1.5));
	fmt.Printf("NextXsquare => %f\n", NextXsquare(3));
	fmt.Printf("NextNegativeBinomial => %d\n", NextNegativeBinomial(.5, 10));
	fmt.Printf("NextStudentsT => %f\n", NextStudentsT(7));
	fmt.Printf("NextF => %f\n", NextF(7, 3));
	fmt.Printf("NextWishart => %v\n",
		NextWishart(100, [][]float64{[]float64{1, 0}, []float64{0, 1}}));
	fmt.Printf("NextInverseWishart => %v\n",
		NextInverseWishart(100, [][]float64{[]float64{1, 0}, []float64{0, 1}}));
}
