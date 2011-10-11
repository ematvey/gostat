package stat

import (
	"math"
	. "gomatrix.googlecode.com/hg/matrix"
)

func MatrixNormal_PDF(M, Omega, Sigma *DenseMatrix) func(A *DenseMatrix) float64 {
	var n float64 = float64(M.Rows())
	var p float64 = float64(Omega.Rows())
	if M.Cols() != int(p) {
		panic("M.Cols != Omega.Rows")
	}
	if Omega.Cols() != int(p) {
		panic("Omega.Cols != Omega.Rows")
	}
	if Sigma.Rows() != int(n) {
		panic("Sigma.Rows != M.Rows")	
	}
	if Sigma.Cols() != int(n) {
		panic("Sigma.Cols != M.Rows")	
	}
	norm := math.Pow(2 * math.Pi, -0.5 * n * p)
	norm *= math.Pow(Omega.Det(), -0.5 * n)
	norm *= math.Pow(Sigma.Det(), -0.5 * p)

	return func(X *DenseMatrix) (p float64) {
		p = norm

		sinv, err := Sigma.Inverse()
		if err != nil {
			panic(err)
		}
		oinv, err := Omega.Inverse()
		if err != nil {
			panic(err)
		}
		diff, err := X.MinusDense(M)
		if err != nil {
			panic(err)
		}
		inner := oinv
		
		inner, err = inner.TimesDense(diff.Transpose())
		if err != nil { panic(err) }

		inner, err = inner.TimesDense(sinv)
		if err != nil { panic(err) }

		inner, err = inner.TimesDense(diff)
		if err != nil { panic(err) }

		innerTrace := inner.Trace()

		p *= math.Exp(-0.5 * innerTrace)

		return
	}
}
func MatrixNormal_LnPDF(M, Omega, Sigma *DenseMatrix) func(A *DenseMatrix) float64 {
	var n float64 = float64(M.Rows())
	var p float64 = float64(Omega.Rows())
	if M.Cols() != int(p) {
		panic("M.Cols != Omega.Rows")
	}
	if Omega.Cols() != int(p) {
		panic("Omega.Cols != Omega.Rows")
	}
	if Sigma.Rows() != int(n) {
		panic("Sigma.Rows != M.Rows")	
	}
	if Sigma.Cols() != int(n) {
		panic("Sigma.Cols != M.Rows")	
	}

	sinv, err := Sigma.Inverse()
	if err != nil {
		panic(err)
	}
	oinv, err := Omega.Inverse()
	if err != nil {
		panic(err)
	}
	
	norm := (2 * math.Pi) * (-0.5 * n * p)
	norm += Omega.Det() * (-0.5 * n)
	norm += Sigma.Det() * (-0.5 * p)

	return func(X *DenseMatrix) (lp float64) {
		lp = norm
		diff, err := X.MinusDense(M)
		if err != nil {
			panic(err)
		}
		inner := oinv
		
		inner, err = inner.TimesDense(diff.Transpose())
		if err != nil { panic(err) }

		inner, err = inner.TimesDense(sinv)
		if err != nil { panic(err) }

		inner, err = inner.TimesDense(diff)
		if err != nil { panic(err) }

		innerTrace := inner.Trace()

		lp += -0.5 * innerTrace

		return
	}
}
func MatrixNormal(M, Omega, Sigma *DenseMatrix) func() (X *DenseMatrix) {
	Mv := Vectorize(M)
	Cov := Kronecker(Omega, Sigma)
	normal := MVNormal(Mv.Array(), Cov.Arrays())
	return func() (X *DenseMatrix) {
		Xv := MakeDenseMatrix(normal(), Mv.Rows(), Mv.Cols())
		X = Unvectorize(Xv, M.Rows(), M.Cols())
		return
	}
}
func NextMatrixNormal(M, Omega, Sigma *DenseMatrix) (X *DenseMatrix) {
	return MatrixNormal(M, Omega, Sigma)()
}

func MVNormal_PDF(m []float64, S [][]float64) func(x []float64) float64 {
	p := len(m)
	μ := MakeDenseMatrix(m, p, 1)
	backμ := μ.DenseMatrix()
	backμ.Scale(-1)
	Σ := MakeDenseMatrixStacked(S)

	Σdet := Σ.Det()
	ΣdetRt := sqrt(Σdet)
	Σinv, _ := Σ.Inverse()

	normalization := pow(2*π, -float64(p)/2) / ΣdetRt

	return func(x []float64) float64 {
		δ, _ := MakeDenseMatrix(x, p, 1).PlusDense(backμ)
		tmp := δ.Transpose()
		tmp, _ = tmp.TimesDense(Σinv)
		tmp, _ = tmp.TimesDense(δ)
		f := tmp.Get(0, 0)
		return normalization * exp(-f/2)
	}
}
func NextMVNormal(μ []float64, Σ [][]float64) []float64 {
	n := len(μ)
	x := Zeros(n, 1)
	for i := 0; i < n; i++ {
		x.Set(i, 0, NextNormal(0, 1))
	}
	C, err := MakeDenseMatrixStacked(Σ).Cholesky()
	Cx, err := C.TimesDense(x)
	μCx, err := MakeDenseMatrix(μ, n, 1).PlusDense(Cx)
	if err != nil {
		panic(err)
	}
	return μCx.Array()
}

func MVNormal(μ []float64, Σ [][]float64) func() []float64 {
	C, _ := MakeDenseMatrixStacked(Σ).Cholesky()
	n := len(μ)
	M := MakeDenseMatrix(μ, n, 1)
	return func() []float64 {
		x := Zeros(n, 1)
		for i := 0; i < n; i++ {
			x.Set(i, 0, NextNormal(0, 1))
		}
		Cx, _ := C.TimesDense(x)
		MCx, _ := M.PlusDense(Cx)
		return MCx.Array()
	}
}

func Wishart_PDF(n int, V [][]float64) func(W [][]float64) float64 {
	Vm := MakeDenseMatrixStacked(V)
	p := Vm.Rows()
	Vdet := Vm.Det()
	Vinv, _ := Vm.Inverse()
	normalization := pow(2, -0.5*float64(n*p)) *
		pow(Vdet, -0.5*float64(n)) /
		Γ(0.5*float64(n))
	return func(W [][]float64) float64 {
		Wm := MakeDenseMatrixStacked(W)
		VinvWm, _ := Vinv.Times(Wm)
		return normalization * pow(Wm.Det(), 0.5*float64(n-p-1)) *
			exp(-0.5*VinvWm.Trace())
	}
}
func Wishart_LnPDF(n int, V [][]float64) func(W [][]float64) float64 {
	Vm := MakeDenseMatrixStacked(V)
	p := Vm.Rows()
	Vdet := Vm.Det()
	Vinv, _ := Vm.Inverse()
	normalization := log(2)*(-0.5*float64(n*p)) +
		log(Vdet)*(-0.5*float64(n)) -
		LnΓ(0.5*float64(n))
	return func(W [][]float64) float64 {
		Wm := MakeDenseMatrixStacked(W)
		VinvWm, _ := Vinv.Times(Wm)
		return normalization +
			log(Wm.Det())*0.5*float64(n-p-1) -
			0.5*VinvWm.Trace()
	}
}
func NextWishart(n int, V [][]float64) [][]float64 {
	return Wishart(n, V)()
}
func Wishart(n int, V [][]float64) func() [][]float64 {
	p := len(V)
	zeros := make([]float64, p)
	rowGen := MVNormal(zeros, V)
	return func() [][]float64 {
		x := make([][]float64, n)
		for i := 0; i < n; i++ {
			x[i] = rowGen()
		}
		X := MakeDenseMatrixStacked(x)
		S, _ := X.Transpose().TimesDense(X)
		return S.Arrays()
	}
}

func InverseWishart_PDF(n int, V [][]float64) func(B [][]float64) float64 {
	p := len(V)
	Ψ := MakeDenseMatrixStacked(V)
	Ψdet := Ψ.Det()
	normalization := pow(Ψdet, -0.5*float64(n)) *
		pow(2, -0.5*float64(n*p)) /
		Γ(float64(n)/2)
	return func(B [][]float64) float64 {
		Bm := MakeDenseMatrixStacked(B)
		Bdet := Bm.Det()
		Binv, _ := Bm.Inverse()
		ΨBinv, _ := Ψ.Times(Binv)
		return normalization *
			pow(Bdet, -.5*float64(n+p+1)) *
			exp(-0.5*ΨBinv.Trace())
	}
}
func InverseWishart_LnPDF(n int, V [][]float64) func(W [][]float64) float64 {
	p := len(V)
	Ψ := MakeDenseMatrixStacked(V)
	Ψdet := Ψ.Det()
	normalization := log(Ψdet)*-0.5*float64(n) +
		log(2)*-0.5*float64(n*p) -
		LnΓ(float64(n)/2)
	return func(B [][]float64) float64 {
		Bm := MakeDenseMatrixStacked(B)
		Bdet := Bm.Det()
		Binv, _ := Bm.Inverse()
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
	p := len(V)
	zeros := make([]float64, p)
	rowGen := MVNormal(zeros, V)
	return func() [][]float64 {
		x := make([][]float64, n)
		for i := 0; i < n; i++ {
			x[i] = rowGen()
		}
		X := MakeDenseMatrixStacked(x)
		S, _ := X.Transpose().TimesDense(X)
		Sinv, _ := S.Inverse()
		return Sinv.Arrays()
	}
}

func Dirichlet_PDF(α []float64) func(θ []float64) float64 {
	return func(θ []float64) float64 {
		if len(θ) != len(α) {
			return 0
		}
		l := float64(1.0)
		totalα := float64(0)
		for i := 0; i < len(θ); i++ {
			if θ[i] < 0 || θ[i] > 1 {
				return 0
			}
			l *= pow(θ[i], α[i]-1)
			l /= Γ(α[i])
			totalα += α[i]
		}
		l *= Γ(totalα)
		return l
	}
}
func Dirichlet_LnPDF(α []float64) func(x []float64) float64 {
	return func(x []float64) float64 {
		if len(x) != len(α) {
			return negInf
		}
		l := fZero
		totalα := float64(0)
		for i := 0; i < len(x); i++ {
			if x[i] < 0 || x[i] > 1 {
				return negInf
			}
			l += (α[i] - 1) * log(x[i])
			l -= LnΓ(α[i])
			totalα += α[i]
		}
		l += LnΓ(totalα)
		return l
	}
}
func NextDirichlet(α []float64) []float64 {
	x := make([]float64, len(α))
	sum := fZero
	for i := 0; i < len(α); i++ {
		x[i] = NextGamma(α[i], 1.0)
		sum += x[i]
	}
	for i := 0; i < len(α); i++ {
		x[i] /= sum
	}
	return x
}
func Dirichlet(α []float64) func() []float64 {
	return func() []float64 { return NextDirichlet(α) }
}
