package stat

import "math"

//Fact(n) = n*Fact(n-1)
func Fact(n int64) int64 {
	return PartialFact(n, 0);
}
//LnFact(n) = log(n)+LnFact(n-1)
func LnFact(n int64) float64 {
	return LnPartialFact(n, 0);
}

//returns Fact(n)/Fact(m)
func PartialFact(n int64, m int64) int64 {
	if n == m {
		return 1;
	}
	return n*PartialFact(n-1, m);
}

//returns LnFact(n)-LnFact(m)
func LnPartialFact(n int64, m int64) float64 {
	if n == m {
		return 0;
	}
	return log(float64(n))+LnPartialFact(n-1, m);
}

func Choose(n int64, i int64) int64 {
	smaller := i;
	if n-i < smaller {
		smaller = n-i;
	}
	return PartialFact(n, smaller)/Fact(smaller);
}

func LnChoose(n int64, i int64) float64 {
	smaller := i;
	if n-i < smaller {
		smaller = n-i;
	}
	return LnPartialFact(n, smaller)-LnFact(smaller);
}

func ChooseMany(i []int64) int64 {
	return 0
}

var lanczos_coef []float64 = []float64{
     0.99999999999980993,
   676.5203681218851,
 -1259.1392167224028,
   771.32342877765313,
  -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
     9.9843695780195716e-6,
     1.5056327351493116e-7 }

//The Gamma function
var Γ = math.Gamma

var sqrt2pi = math.Sqrt(2 * math.Pi)

func LnΓ(x float64) float64 {	
	tmp := (x - 0.5) * log(x + 4.5) - (x + 4.5)
	ser := 1.0 +
			76.1800917300 / (x + 0) - 86.5053203300 / (x + 1) +
			24.0140982200 / (x + 2) - 1.23173951600 / (x + 3) +
			0.00120858003 / (x + 4) - 0.00000536382 / (x + 5)
	return tmp + log(ser * sqrt2pi)
}

func B(x float64, y float64) float64 {
	return Γ(x)*Γ(y)/Γ(x+y)
}
func LnB(x float64, y float64) float64 {
	return LnΓ(x)+LnΓ(y)-LnΓ(x+y)
}

/*

	public static double logGammaP(int p, double x) {
		double r = p*(p-1)*.25*Math.log(Math.PI);
		
		for (int j=1; j<=p; j++) {
			r += logGamma(x+.5*(1-j));
		}
		
		return r;
	}
	public static double logGammaPRatio(int p, double numerator, double denominator) {
		double r = 0;

		for (int j=1; j<=p; j++) {
			r += logGamma(numerator+.5*(1-j));
			r -= logGamma(denominator+.5*(1-j));
		}
		
		return r;
	}
*/

func LnΓp(p int, x float64) (r float64) {
	pf := float64(p)
	r = pf*(pf-1)*.25*math.Log(math.Pi)
	for j:=float64(1); j<=pf; j++ {
		r += LnΓ(x+.5*(1-j))
	}
	return
}

//LnΓp(x)/LnΓp(y)
func LnΓpRatio(p int, x, y float64) (r float64) {
	pf := float64(p)
	for j:=float64(1); j<=pf; j++ {
		r += LnΓ(x+.5*(1-j))
		r -= LnΓ(y+.5*(1-j))
	}
	return
}

