package stat

import (
	. "gomatrix.googlecode.com/hg/matrix";
)

func MVNormal_PDF(m []float64, S [][]float64) func(x []float64) float64 {
	p := len(m);
	μ := MakeDenseMatrix(m, p, 1);
	backμ := μ.DenseMatrix()
	backμ.Scale(-1)
	Σ := MakeDenseMatrixStacked(S);
	
	Σdet := Σ.Det();
	ΣdetRt := sqrt(Σdet);
	Σinv, _ := Σ.Inverse();
	
	normalization := pow(2*π, -float64(p)/2)/ΣdetRt;
	
	return func(x []float64) float64 {
		δ, _ := MakeDenseMatrix(x, p, 1).PlusDense(backμ);
		tmp := δ.Transpose()
		tmp, _ = tmp.TimesDense(Σinv)
		tmp, _ = tmp.TimesDense(δ)
		f := tmp.Get(0, 0);
		return normalization*exp(-f/2);
	}
}
func NextMVNormal(μ []float64, Σ [][]float64) []float64 {
	n := len(μ);
	x := Zeros(n, 1);
	for i:=0; i<n; i++ {
		x.Set(i, 0, NextNormal(0, 1));
	}
	C, _ := MakeDenseMatrixStacked(Σ).Cholesky();
	Cx, _ := C.TimesDense(x)
	μCx, _ := MakeDenseMatrix(μ, n, 1).PlusDense(Cx)
	return μCx.Array()
}

func MVNormal(μ []float64, Σ [][]float64) func() []float64 {
	C, _ := MakeDenseMatrixStacked(Σ).Cholesky();
	n := len(μ);
 	M := MakeDenseMatrix(μ, n, 1);
	return func() []float64 {
		x := Zeros(n, 1);
		for i:=0; i<n; i++ {
			x.Set(i, 0, NextNormal(0, 1));
		}
		Cx, _ := C.TimesDense(x)
		MCx, _ := M.PlusDense(Cx)
		return MCx.Array()
	}
}

func Wishart_PDF(n int, V [][]float64) func(W [][]float64) float64 {
	Vm := MakeDenseMatrixStacked(V);
	p := Vm.Rows();
	Vdet := Vm.Det();
	Vinv, _ := Vm.Inverse();
	normalization := pow(2, -0.5*float64(n*p))*
					pow(Vdet, -0.5*float64(n))/
					Γ(0.5*float64(n));
	return func(W [][]float64) float64 {
		Wm := MakeDenseMatrixStacked(W);
		VinvWm, _ := Vinv.Times(Wm)
		return normalization*pow(Wm.Det(), 0.5*float64(n-p-1)) *
			exp(-0.5*VinvWm.Trace());
	}
}
func Wishart_LnPDF(n int, V [][]float64) func(W [][]float64) float64 {
	Vm := MakeDenseMatrixStacked(V);
	p := Vm.Rows();
	Vdet := Vm.Det();
	Vinv, _ := Vm.Inverse();
	normalization := log(2)*(-0.5*float64(n*p)) +
		log(Vdet)*(-0.5*float64(n)) -
		LnΓ(0.5*float64(n));
	return func(W [][]float64) float64 {
		Wm := MakeDenseMatrixStacked(W);
		VinvWm, _ := Vinv.Times(Wm)
		return normalization +
			log(Wm.Det())*0.5*float64(n-p-1) -
			0.5*VinvWm.Trace();
	}
}
func NextWishart(n int, V [][]float64) [][]float64 {
	return Wishart(n, V)()
}
func Wishart(n int, V [][]float64) func() [][]float64 {
	p := len(V);
	zeros := make([]float64, p);
	rowGen := MVNormal(zeros, V);
	return func() [][]float64 {			
		x := make([][]float64, n);
		for i:=0; i<n; i++ {
			x[i] = rowGen()
		}
		X := MakeDenseMatrixStacked(x);
		S, _ := X.Transpose().TimesDense(X);
		return S.Arrays();
	}
}

func InverseWishart_PDF(n int, V [][]float64) func(B [][]float64) float64 {
	p := len(V);
	Ψ := MakeDenseMatrixStacked(V);
	Ψdet := Ψ.Det();
	normalization := pow(Ψdet, -0.5*float64(n)) *
		pow(2, -0.5*float64(n*p)) /
		Γ(float64(n)/2);
	return func(B [][]float64) float64 {
		Bm := MakeDenseMatrixStacked(B);
		Bdet := Bm.Det();
		Binv, _ := Bm.Inverse();
		ΨBinv, _ := Ψ.Times(Binv)
		return normalization *
			pow(Bdet, -.5*float64(n+p+1)) *
			exp(-0.5*ΨBinv.Trace())
	}
}
func InverseWishart_LnPDF(n int, V [][]float64) func(W [][]float64) float64 {
	p := len(V);
	Ψ := MakeDenseMatrixStacked(V);
	Ψdet := Ψ.Det();
	normalization := log(Ψdet)*-0.5*float64(n) +
		log(2)*-0.5*float64(n*p) -
		LnΓ(float64(n)/2);
	return func(B [][]float64) float64 {
		Bm := MakeDenseMatrixStacked(B);
		Bdet := Bm.Det();
		Binv, _ := Bm.Inverse();
		ΨBinv, _ := Ψ.Times(Binv)
		return normalization +
			log(Bdet)*-.5*float64(n+p+1) +
			-0.5*ΨBinv.Trace()
	}
}
func NextInverseWishart(n int, V [][]float64) [][]float64 {
	return InverseWishart(n, V)()
}
func InverseWishart(n int, V [][]float64) func() [][]float64 {
	p := len(V);
	zeros := make([]float64, p);
	rowGen := MVNormal(zeros, V);
	return func() [][]float64 {
		x := make([][]float64, n);
		for i:=0; i<n; i++ {
			x[i] = rowGen();
		}
		X := MakeDenseMatrixStacked(x);
		S, _ := X.Transpose().TimesDense(X);
		Sinv, _ := S.Inverse()
		return Sinv.Arrays();
	}
}

func Dirichlet_PDF(α []float64) func(θ []float64) float64 {
	return func(θ []float64) float64 {
		if len(θ) != len(α) {
			return 0
		}
		l := float64(1.0);
		totalα := float64(0);
		for i := 0; i < len(θ); i++ {
			if θ[i] < 0 || θ[i] > 1 {
				return 0
			}
			l *= pow(θ[i], α[i]-1);
			l *= Γ(α[i]);
		}
		l /= Γ(totalα);
		return l;
	}
}
func Dirichlet_LnPDF(α []float64) func(x []float64) float64 {
	return func(x []float64) float64 {
		if len(x) != len(α) {
			return negInf
		}
		l := fZero;
		totalα := float64(0.0);
		for i := 0; i < len(x); i++ {
			if x[i] < 0 || x[i] > 1 {
				return negInf
			}
			l += (α[i] - 1) * log(x[i]);
			l += LnΓ(α[i]);
		}
		l -= LnΓ(totalα);
		return l;
	}
}
func NextDirichlet(α []float64) []float64 {
	x := make([]float64, len(α));
	sum := fZero;
	for i := 0; i < len(α); i++ {
		x[i] = NextGamma(α[i], 1.0);
		sum += x[i];
	}
	for i := 0; i < len(α); i++ {
		x[i] /= sum
	}
	return x;
}
func Dirichlet(α []float64) func() []float64 {
	return func() []float64 { return NextDirichlet(α) }
}
