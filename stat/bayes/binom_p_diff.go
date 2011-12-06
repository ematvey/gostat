/*
Bayesian inference for the difference between two binomial proportions using the Normal (Gaussian) approximation.
Bolstad 2007 (2e): Chapter 13.4: 248-249.

Compare the proportions of certain attribute in two populations. the true proportions are pi1 and pi2, unknown.
We take a random sample from each of the populations and observe y1, y2 ... number of ibnstances having the attribute.
The distribution y1|pi1 is binomial(n1, pi1), similarly for y2|pi2, and they are independent.
Let there be independent priors for pi1 ... beta(a1, b1), and similarly for pi2.
Posterior for pi1 is beta(a1_post, b1_post), where a1_post = a1 +y1, and b1_post = b1 + n1 -y1, similarly for pi2.
Approximate each posterior distribution with normal distribution having the same mean and variance as the beta.
The posterior of pi_d = pi1 - pi2 is approximately normal(md_post, vard_post), where:
md_post = a1_post/(a1_post+b1_post) - a2_post/(a2_post+b2_post), and
vard_post = a1_post*b1_post/SQR(a1_post+b1_post)*(a1_post+b1_post+1)  +  a2_post*b2_post/SQR(a2_post+b2_post)*(a2_post+b2_post+1) 
*/

package bayes
import (
	"fmt"
	"math"
	s "gostat.googlecode.com/hg/stat"
)


/* 
Mean of posterior distribution of unknown difference of binomial proportions, approximated by Normal distribution
Bolstad 2007 (2e): 248.
untested ...
*/
func BinomDiffPropNormApproxMean(a1, b1, a2, b2 float64, n1, n2 int64) float64 {

	a1_post := a1 +y1
	b1_post := b1 + n1 -y1
	a2_post := a2 +y2
	b2_post := b2 + n2 -y2

	md_post := a1_post/(a1_post+b1_post) - a2_post/(a2_post+b2_post)
	return(md_post)
}

/* 
Variance of posterior distribution of unknown difference of binomial proportions, approximated by Normal distribution
Bolstad 2007 (2e): 248.
untested ...
*/
func BinomDiffPropNormApproxVar(a1, b1, a2, b2 float64, n1, n2 int64) float64 {

	a1_post = a1 +y1
	b1_post = b1 + n1 -y1
	a2_post = a2 +y2
	b2_post = b2 + n2 -y2

	vard_post = a1_post*b1_post/SQR(a1_post+b1_post)*(a1_post+b1_post+1)  +  a2_post*b2_post/SQR(a2_post+b2_post)*(a2_post+b2_post+1)
	return(vard_post)
}

/*
Credible interval for difference between binomial proportions, approximated by Normal distribution
Bolstad 2007 (2e): 248, eq. 13.13
post_diff_mu = binom_diff_prop_norm_approx_mu()
post_diff_sigma = sqrt(binom_diff_prop_norm_approx_var())
untested ...
*/
func BinomDiffPropCrI(post_diff_mu, post_diff_sigma, alpha float64) (float64, float64) {
	/*
	post_diff_mu		posterior mean for difference of normal means
	post_diff_sigma	posterior standard deviation for difference of normal means
	alpha			posterior probability that the true mean lies outside the credible interval
	*/

	var z float64

	z = s.ZInv_CDF_For(alpha/2)

	low = post_diff_mu - z*post_diff_sigma
	high = post_diff_mu + z*post_diff_sigma
	return low, high
}

/*
One-sided test for difference between binomial proportions, approximated by Normal distribution
Bolstad 2007 (2e): 248-249, eq. 13.14
H0: mu_d <= 0 vs H1: mu_d > 0
Note: The alternative is in the direction we wish to detect, and is what we want to detect.
*/
func BinomDiffPropOneSidedProb(post_diff_mu, post_diff_sigma float64) float64 {
	var prob float64
//	prob =gsl_cdf_ugaussian_P(-post_diff_mu / post_diff_sigma)
	prob = s.Z_CDF_For(-post_diff_mu / post_diff_sigma)
	return(prob)
}

/*
Two-sided test for difference between binomial proportions, approximated by Normal distribution  ///// check it vs the book!!!
Bolstad 2007 (2e): 249
H0: mu_1 - mu_2 == 0 vs H1: mu_1 - mu_2 != 0
func BinomDiffPropTwoSidedProb(post_diff_mu, post_diff_sigma, alpha){
	low, high = norm_diff_means_know_CrI(post_diff_mu, post_diff_sigma, alpha)

	if 0 < low || 0 > high return(REJECT) else return(ACCEPT)
}
*/

