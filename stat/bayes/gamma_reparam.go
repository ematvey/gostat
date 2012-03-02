// Parameterizing a Gamma distribution by mode and sd.
// It is more intuitive to start with the mode and standard deviation, instead of the mean and standard deviation as used in the Kruschke (2011) book. 
// The reason is that the gamma distribution is typically very skewed, and therefore the location of the mean is not very intuitive. 
// This function computes the shape and rate parameters of the gamma distribution from a desired mode and standard deviation.
// After http://doingbayesiandataanalysis.blogspot.com/2012/01/parameterizing-gamma-distribution-by.html
package bayes

import (
	"math"
)

// Parameterizing a Gamma distribution by mode and standard deviation
func GammaMSParam(mode, sd float64) (rate, shape float64) {
	rate = (mode + math.Sqrt(mode*mode+4*sd*sd)) / (2 * sd * sd)
	shape = 1 + mode*rate
	return
}
