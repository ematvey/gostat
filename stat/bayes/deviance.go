// Bayesian deviance
// Aitkin 2010: 41

package bayes
import "math"

// Bayesian deviance
func deviance(likelihood float64) float64 {
	return -2*math.Log(likelihood)
}

// Bayesian deviance difference
func dev_diff(like1, like2 float64) float64 {
	lr := like1/like2
	return -2*math.Log(lr)
}


