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
func Γ(z float64) float64 {
	const g = 7;
	if z < 0.5 {
		return π/(math.Sin(π*z)*Γ(1-z));
	}
	z -= 1;
	x := lanczos_coef[0];
	for i:=1; i<g+2; i++ {
		x += lanczos_coef[i]/(z+float64(i));
	}
	t := z+g+.5;
	return sqrt(2*π)*pow(t, z+.5)*exp(-1*t)*x;
}

func LnΓ(x float64) float64 {
	f := fZero;
	var z float64;
	
	if (x<7) {
		f=1;
		z=x-1;
		for ; z<7; z++ {
			f*=z;
		}
		x=z;
		f=-log(f);
	}
	z = 1/(x*x);
	return f + (x-0.5)*log(x) - x + .918938533204673 +
			(((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z +
			.083333333333333)/x;  
}

func B(x float64, y float64) float64 {
	return Γ(x)*Γ(y)/Γ(x+y)
}
func LnB(x float64, y float64) float64 {
	return LnΓ(x)+LnΓ(y)-LnΓ(x+y)
}
