// Bayesian inference about the parameter vector p of multinomial distribution.
// Conjugate prior is Dirichlet(|α|), conjugate posterior is Dirichlet(|α+x|).

package bayes
import (
	"fmt"
	. "gostat.googlecode.com/hg/stat"
)

// Posterior PDF, Dirichlet prior
func MultinomPi_PDF_DirPri(α, x []float64) float64 {
	// if α == nil, use ones
	if α == nil {
		for i := 0; i < len(x); i++ {
			α[i] = 1
		}
	}

	if len(α) != len(x){
		panic(fmt.Sprintf("len(α) != len(x)"))
	}

	for i := 0; i < len(x); i++ {
		α[i] =  α[i]+x[i]	// posterior params
	}
	return	Dirichlet_PDF_At(α, x)
}

// Sampling from posterior, Dirichlet prior
// Returns an array of sampled Multinomial Pi's
func NextMultinomPi(α, x []float64) []float64 {
	for i := 0; i < len(x); i++ {
		α[i] =  α[i]+(x[i])	// posterior params
	}
	return NextDirichlet(α)
}

