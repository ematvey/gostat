package stat

import (
	"math"
	"gomatrix.googlecode.com/hg/matrix"
)

func MatrixT_PDF(M, Sigma, Omega *matrix.DenseMatrix, n int) func(T *matrix.DenseMatrix) (l float64) {
	nf := float64(n)
	p := M.Rows()
	pf := float64(p)
	m := M.Cols()
	mf := float64(m)

	if m != Sigma.Rows() || m != Sigma.Cols() {
		panic("Sigma is not m x m")
	}
	if p != Omega.Rows() || p != Omega.Cols() {
		panic("Omega is not p x p")
	}
	if n <= 0 {
		panic("n <= 0")
	}

	var norm float64 = 1

	norm *= GammaPRatio(p, 0.5 * (nf + mf + pf - 1), 0.5 * (nf + pf - 1))
	norm *= math.Pow(math.Pi, -0.5 * mf * pf)
	norm *= math.Pow(Sigma.Det(), -0.5 * mf)
	norm *= math.Pow(Omega.Det(), -0.5 * pf)
	
	SigmaInv, err := Sigma.Inverse()
	if err != nil { panic(err) }
	OmegaInv, err := Omega.Inverse()
	if err != nil { panic(err) }

	return func(T *matrix.DenseMatrix) (l float64) {
		l = norm

		diff, err := T.MinusDense(M)
		if err != nil { panic(err) }
		inner := SigmaInv.Copy()
		inner, _ = inner.TimesDense(diff)
		inner, _ = inner.TimesDense(OmegaInv)
		inner, _ = inner.TimesDense(diff.Transpose())

		l *= math.Pow(inner.Det(), -0.5 * (nf + mf + pf - 1))

		return
	}
}

func MatrixT_LnPDF(M, Sigma, Omega *matrix.DenseMatrix, n int) func(T *matrix.DenseMatrix) (ll float64) {
	nf := float64(n)
	p := M.Rows()
	pf := float64(p)
	m := M.Cols()
	mf := float64(m)

	if m != Sigma.Rows() || m != Sigma.Cols() {
		panic("Sigma is not m x m")
	}
	if p != Omega.Rows() || p != Omega.Cols() {
		panic("Omega is not p x p")
	}
	if n <= 0 {
		panic("n <= 0")
	}

	var norm float64 = 0

	norm += LnGammaPRatio(p, 0.5 * (nf + mf + pf - 1), 0.5 * (nf + pf - 1))
	norm += math.Pow(math.Pi, -0.5 * mf * pf)
	norm += math.Pow(Sigma.Det(), -0.5 * mf)
	norm += math.Pow(Omega.Det(), -0.5 * pf)
	
	SigmaInv, err := Sigma.Inverse()
	if err != nil { panic(err) }
	OmegaInv, err := Omega.Inverse()
	if err != nil { panic(err) }

	return func(T *matrix.DenseMatrix) (ll float64) {
		ll = norm

		diff, err := T.MinusDense(M)
		if err != nil { panic(err) }
		inner := SigmaInv.Copy()
		inner, _ = inner.TimesDense(diff)
		inner, _ = inner.TimesDense(OmegaInv)
		inner, _ = inner.TimesDense(diff.Transpose())

		ll += math.Log(inner.Det()) * -0.5 * (nf + mf + pf - 1)

		return
	}
}

func MatrixT(M, Sigma, Omega *matrix.DenseMatrix, n int) func() (T *matrix.DenseMatrix) {
	p := M.Rows()
	m := M.Cols()

	if m != Sigma.Rows() || m != Sigma.Cols() {
		panic("Sigma is not m x m")
	}
	if p != Omega.Rows() || p != Omega.Cols() {
		panic("Omega is not p x p")
	}
	if n <= 0 {
		panic("n <= 0")
	}

	SigmaInv, err := Sigma.Inverse()
	if err != nil { panic(err) }

	Sdist := Wishart(n+p-1, SigmaInv)

	Xdist := MatrixNormal(matrix.Zeros(p, p), matrix.Eye(p), Omega)

	return func() (T *matrix.DenseMatrix) {
		S := Sdist()
		Sc, err := S.Cholesky()
		if err != nil { panic(err) }
		X := Xdist()
		T, err = Sc.Transpose().TimesDense(X)
		err = T.AddDense(M)
		if err != nil { panic(err) }
		return
	}
}

func NextMatrixT(M, Sigma, Omega *matrix.DenseMatrix, n int) (T *matrix.DenseMatrix) {
	return MatrixT(M, Sigma, Omega, n)()
}
