package bayes

import (
	"math"
)

func dsign(a, b float64) float64 {
	if b < 0 { // if negative
		a = -math.Abs(a)
	} else {
		a = math.Abs(a)
	}
	return a
}

// just for testing purposes
func parab(x float64) float64 {
	a := 1.0 / 3.0
	return (x - a) * (x - a)
}

// Fmin returns an approximation  x  to the point where  f  attains a minimum  on
// the interval  (ax,bx).
// params:
// ax    left endpoint of initial interval
// bx    right endpoint of initial interval
// f     function subprogram which evaluates  f(x)  for any  x
//	 in the interval  (ax,bx)
// tol   desired length of the interval of uncertainty of the final
//	 result
// output:
// abcissa approximating the point where  f  attains a minimum.
// The method used is a combination of  golden  section  search  and
// successive parabolic interpolation.  convergence is never much slower
// than  that  for  a  fibonacci search.  if  f  has a continuous second
// derivative which is positive at the minimum (which is not  at  ax  or
// bx),  then  convergence  is  superlinear, and usually of the order of
// about  1.324....
//     the function  f  is never evaluated at two points closer together
// than  eps*abs(fmin) + (tol/3), where eps is  approximately the square
// root  of  the  relative  machine  precision.   if   f   is a unimodal
// function and the computed values of   f   are  always  unimodal  when
// separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
// the abcissa of the global minimum of  f  on the interval  ax,bx  with
// an error less than  3*eps*abs(fmin) + tol.  if   f   is not unimodal,
// then fmin may approximate a local, but perhaps non-global, minimum to
// the same accuracy.
// This function is a slightly modified  version  of  the
// Algol60 procedure  'localmin'  given in Richard Brent, Algorithms for
// minimization without derivatives, Prentice - Hall (1973).
// Based on http://www.netlib.org/fmm/fmin.f
func fmin(f func(float64) float64, ax, bx, tol float64) float64 {

	var (
		a, b, c, d, e, eps, xm, p, q, r, tol1, tol2, u, v, w, fu, fv, fw, fx, x float64
	)

	gsNeeded := false
	// c is the squared inverse of the golden ratio
	c = 0.5 * (3. - math.Sqrt(5))

	// eps is approximately the square root of the relative machine precision
	eps = 1.0

	for {
		eps = eps / 2.0
		tol1 = 1.0 + eps
		if tol1 <= 1.00 {
			break
		}
	}

	eps = math.Sqrt(eps)
	a = ax
	b = bx
	v = a + c*(b-a)
	w = v
	x = v
	e = 0.0
	fx = f(x)
	fv = fx
	fw = fx

	for {
		xm = 0.5 * (a + b)
		tol1 = eps*math.Abs(x) + tol/3.0
		tol2 = 2.0 * tol1

		// check stopping criterion
		if math.Abs(x-xm) <= (tol2 - 0.5*(b-a)) {
			break
		}

		// is golden-section necessary?
		if math.Abs(e) > tol1 {

			// fit parabola
			r = (x - w) * (fx - fv)
			q = (x - v) * (fx - fw)
			p = (x-v)*q - (x-w)*r
			q = 2.00 * (q - r)
			if q > 0.0 {
				p = -p
			}
			q = math.Abs(q)
			r = e
			e = d

			// is parabola acceptable? If not, do a golden section
			if math.Abs(p) < math.Abs(0.5*q*r) && p > q*(a-x) && p < q*(b-x) { // go to 40

				// a parabolic interpolation step
				d = p / q
				u = x + d
				// f must not be evaluated too close to ax or bx
				if (u - a) < tol2 {
					d = dsign(tol1, xm-x)
				}
				if (b - u) < tol2 {
					d = dsign(tol1, xm-x)
				}
				//	go to 50
			} else {
				gsNeeded = true
			}

		} else {
			gsNeeded = true
		}
		if gsNeeded {
			// do a golden-section

			if x >= xm {
				e = a - x
			} else {
				e = b - x
			}
			d = c * e
		}

		// f must not be evaluated too close to x
		//   50 
		if math.Abs(d) >= tol1 {
			u = x + d
		}
		if math.Abs(d) < tol1 {
			u = x + dsign(tol1, d)
		}
		fu = f(u)

		// update  a, b, v, w, and x
		if !(fu > fx) { //go to 60
			if u >= x {
				a = x
			}
			if u < x {
				b = x
			}
			v = w
			fv = fw
			w = x
			fw = fx
			x = u
			fx = fu
			continue
		} else {
			if u < x {
				a = u
			} else {
				b = u
			}
		}
		if (fu <= fw) || (w == x) { //go to 70
			v = w
			fv = fw
			w = u
			fw = fu
		} else if fu <= fv || v == x || v == w { // go to 80
			v = u
			fv = fu
		} // if
	} // for
	return x
}
