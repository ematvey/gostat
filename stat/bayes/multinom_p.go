// Bayesian inference about the parameter vector p of multinomial distribution.
// Conjugate prior is Dirichlet(|α|), conjugate posterior is Dirichlet(|α+x|).
// http://en.wikipedia.org/wiki/Dirichlet_distribution#Conjugate_to_multinomial

package bayes

import (
	"fmt"
	. "gostat.googlecode.com/hg/stat/prob"
)

// Posterior PDF, Dirichlet prior
// for Haldane improper prior, use α[i] = 0
// Ericson 1969 recommends prior with sum(α[i]) small, of the order of 1
// Aitkin 2010:96-107
func MultinomPi_PDF_DirPri(α, x []float64) float64 {
	// if α == nil, use Haldane
	if α == nil {
		for i := 0; i < len(x); i++ {
			α[i] = 0
		}
	}

	if len(α) != len(x) {
		panic(fmt.Sprintf("len(α) != len(x)"))
	}

	for i := 0; i < len(x); i++ {
		α[i] += x[i] // posterior params
	}
	return Dirichlet_PDF_At(α, x)
}

// Sampling from posterior, Dirichlet prior
// Returns an array of sampled Multinomial Pi's
func NextMultinomPi(α, x []float64) []float64 {
	for i := 0; i < len(x); i++ {
		α[i] += x[i] // posterior params
	}
	return NextDirichlet(α)
}
