// Student's t distribution
package prob

import (
	. "go-fn.googlecode.com/hg/fn"
	"math"
)

// Student's t distribution: probability density function
func StudentsT_PDF(ν float64) func(x float64) float64 {
	normalization := Γ((ν+1)/2) / (sqrt(ν*π) * Γ(ν/2))
	return func(x float64) float64 {
		return normalization * pow(1+x*x/ν, -(ν+1)/2)
	}
}

// Student's t distribution: logarithm of the probability density function
func StudentsT_LnPDF(ν float64) func(x float64) float64 {
	normalization := LnΓ((ν+1)/2) - log(sqrt(ν*π)) - LnΓ(ν/2)
	return func(x float64) float64 {
		return normalization + log(1+x*x/ν)*-(ν+1)/2
	}
}

// Student's t distribution: random value
// StudentsT(ν) => N(0, 1)*sqrt(ν/NextGamma(ν/2, 2))
func NextStudentsT(ν float64) float64 {
	return NextNormal(0, 1) * sqrt(ν/NextGamma(ν/2, 2))
}

func StudentsT(ν float64) func() float64 {
	return func() float64 {
		return NextStudentsT(ν)
	}
}

// Student's t distribution: cumulative density function
func StudentsT_CDF(ν float64) func(x float64) float64 {
	return func(x float64) float64 {
		var p float64
		if ν <= 0 {
			panic("ν <= 0")
		}

		nx := 1 + (x/ν)*x
		if nx > 1e100 { /* <==>  x*x > 1e100 * ν  */
			/* Danger of underflow. So use Abramowitz & Stegun 26.5.4
			   pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
			   with z = 1/nx,  a = ν/2,  b= 1/2 :
			*/
			p = -0.5*ν*(2*math.Log(math.Abs(x))-math.Log(ν)) - LnB(0.5*ν, 0.5) - math.Log(0.5*ν)
			p = math.Exp(p)
		} else {
			if ν > x*x {
				α := 0.5
				β := ν / 2
				pbeta := Beta_CDF(α, β)
				p = 1 - pbeta(x*x/(ν+x*x))
			} else {
				α := ν / 2
				β := 0.5
				pbeta := Beta_CDF(α, β)
				p = pbeta (1 / nx)
			}
		}

		p /= 2
		if x > 0 {
			p = 1-p
		}
		return p
	}
}

// Student's t distribution: quantile function
// Hill, G.W (1970) "Algorithm 396: Student's t-quantiles"
// CACM 13(10), 619-620.
// Using expm1() takes care of  Lozy (1979) "Remark on Algo.", TOMS
// Applies 2-term Taylor expansion as in Hill, G.W (1981) "Remark on Algo.396", ACM TOMS 7, 250-1
// Improved formula for decision when 1 < df < 2
func StudentsT_Qtl(ν float64) func(p float64) float64 {
	return func(p float64) float64 {
		const eps = 1.e-12
		var q float64
		neg := false
		p_ok := false

		if ν <= 0 || p <0 || p > 1 {
			panic("bad params")
		}

		/*
			    if (ν < 1) { // based on qnt
				const static double accu = 1e-13;
				const static double Eps = 1e-11; // must be > accu

				double ux, lx, nx, pp;

				int iter = 0;

				p = R_DT_qIv(p);

				// Invert pt(.) :
				// 1. finding an upper and lower bound
				if(p > 1 - math.SmallestNonzeroFloat64) return ML_POSINF;
				pp = fmin2(1 - math.SmallestNonzeroFloat64, p// (1 + Eps));
				for(ux = 1.; ux < DBL_MAX && pt(ux, ν, TRUE, FALSE) < pp; ux//= 2);
				pp = p// (1 - Eps);
				for(lx =-1.; lx > -DBL_MAX && pt(lx, ν, TRUE, FALSE) > pp; lx//= 2);

				// 2. interval (lx,ux)  halving
				   regula falsi failed on qt(0.1, 0.1)

				do {
				    nx = 0.5// (lx + ux);
				    if (pt(nx, ν, TRUE, FALSE) > p) ux = nx; else lx = nx;
				} while ((ux - lx) / math.Abs(nx) > accu && ++iter < 1000);

				if(iter >= 1000) ML_ERROR(ME_PRECISION, "qt");

				return 0.5// (lx + ux);
			    }
		*/

		if ν > 1e20 {
			q = Z_Qtl_For(p)
		} else {
			if p < 0.5 {
				neg = true
			}
			if neg {
				p = 2 * p
			} else {
				p = 2 * (0.5 - p + 0.5)
			}

			if math.Abs(ν-2) < eps { // df ~= 2
				if p > math.SmallestNonzeroFloat64 {
					if 3*p < math.SmallestNonzeroFloat64 { // p ~= 0
						q = 1 / math.Sqrt(p)
					} else if p > 0.9 { // p ~= 1
						q = (1 - p) * math.Sqrt(2/(p*(2-p)))
					} else { // eps/3 <= p <= 0.9
						q = math.Sqrt(2/(p*(2-p)) - 2)
					}
				} else { // p << 1, q = 1/math.Sqrt(p) = ...
					panic("q = +Inf")

				}
			} else if ν < 1+eps { // df ~= 1  (df < 1 excluded above): Cauchy
				if p > 0 {
					q = 1 / math.Tan(p*π/2) // == - math.Tan((p+1) * π/2) -- suffers for p ~= 0

				} else { // p = 0, but maybe = 2*exp(p) !
					panic("q = +Inf")
				}
			} else { //-- usual case;  including, e.g.,  df = 1.1
				x := 0.0
				y := 0.0
				a := 1 / (ν - 0.5)
				b := 48 / (a * a)
				c := ((20700*a/b-98)*a-16)*a + 96.36
				d := ((94.5/(b+c)-3)/b + 1) * math.Sqrt(a*π/2) * ν

				y = math.Pow(d*p, 2/ν)
				if y >= math.SmallestNonzeroFloat64 {
					p_ok = true
				}
				if (ν < 2.1 && p > 0.5) || y > 0.05+a { // p > p0(df)
					// Asymptotic inverse expansion about normal
					x = Normal_Qtl_For(0, 1, 0.5*p)

					y = x * x
					if ν < 5 {
						c += 0.3 * (ν - 4.5) * (x + 0.6)
					}
					c = (((0.05*d*x-5)*x-7)*x-2)*x + b + c
					y = (((((0.4*y+6.3)*y+36)*y+94.5)/c-y-3)/b + 1) * x
					y = expm1(a * y * y)
					q = math.Sqrt(ν * y)
				} else { // re-use 'y' from above

					if !p_ok && x < -0.5*log(math.SmallestNonzeroFloat64) { // 0.5* log(math.SmallestNonzeroFloat64)
						// y above might have underflown
						q = math.Sqrt(ν) * exp(-x)
					} else {
						y = ((1/(((ν+6)/(ν*y)-0.089*d-0.822)*(ν+2)*3)+0.5/(ν+4))*y-1)*(ν+1)/(ν+2) + 1/y
						q = math.Sqrt(ν * y)
					}
				}

				// Now apply 2-term Taylor expansion improvement (1-term = Newton): as by Hill (1981) [ref.above]
				dt := StudentsT_PDF(ν)
				pt := StudentsT_CDF(ν)
				for it := 0; it < 10; it++ {
					y = dt(q)
					if y <= 0 {
						break
					}
					x = (1 - pt(q) - p/2) / y
					if math.Abs(x) > 1e-14*math.Abs(q) {
						// Newton (=Taylor 1 term):
						//  q += x 
						// Taylor 2-term : 
						q += x * (1. + x*q*(ν+1)/(2*(q*q+ν)))
					}
				}
			}
		}
		if neg {
			q = -q
		}
		return q
	}
}



