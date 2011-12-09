package stat
import "math"

// Inverse Gamma distribution: probability density function
func InvGamma_PDF(a, b float64) func(x float64) float64 {
	return func(x float64) float64 {
		return math.Exp(a*math.Log(b) - LnΓ(a) - (a+1)*math.Log(x) - b*1.0/x)
	}
}

// Inverse Gamma distribution: probability density function at x
func InvGamma_PDF_At(a, b float64) func(x float64) float64 {
	return func(x float64) float64 {
		return math.Exp(a*math.Log(b) - LnΓ(a) - (a+1)*math.Log(x) - b*1.0/x)
	}
}

// Inverse Gamma distribution: cumulative distribution function
// not to be confused with Inverse CDF of Gamma distribution
func InvGamma_CDF(a, b float64) func(x float64) float64 {
	return func(x float64) float64 {
		1 - func IΓ(a, b*1.0/x)
	}
}


