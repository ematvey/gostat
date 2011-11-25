package inf

/*
 functions for Bayesian linear regression
*/

import (
	"fmt"
	mx "gomatrix.googlecode.com/hg/matrix"
	"gostat.googlecode.com/hg/stat"
)

type KnownVarianceLRPosterior struct {
	Sigma, M, Phi *mx.DenseMatrix

	XXt, YXt *mx.DenseMatrix

	sampler func() *mx.DenseMatrix
}

/*
 M is r x c, o x i
 Sigma is r x r, o x o
 Phi is c x c, i x i

 Sigma matches Y o x 1 output dimension
 Phi matches X i x 1 input dimension
*/
func NewKnownVarianceLRPosterior(M, Sigma, Phi *mx.DenseMatrix) (this *KnownVarianceLRPosterior) {
	if M.Rows() != Sigma.Rows() {
		panic("M.Rows != Sigma.Rows")
	}
	if M.Cols() != Phi.Cols() {
		panic("M.Cols != Phi.Cols")
	}
	if Sigma.Rows() != Sigma.Cols() {
		panic("Sigma is not square")
	}
	if Phi.Rows() != Phi.Cols() {
		panic("Phi is not square")
	}
	this = &KnownVarianceLRPosterior {
		M: M,
		Sigma: Sigma,
		Phi: Phi,
		XXt: mx.Zeros(Phi.Cols(), Phi.Cols()),
		YXt: mx.Zeros(Sigma.Cols(), Phi.Cols()),
	}

	return
}

func (this *KnownVarianceLRPosterior) Insert(x, y *mx.DenseMatrix) {
	xxt, _ := x.TimesDense(x.Transpose())
	this.XXt.Add(xxt)
	yxt, _ := y.TimesDense(x.Transpose())
	this.YXt.Add(yxt)
}

func (this *KnownVarianceLRPosterior) Remove(x, y *mx.DenseMatrix) {
	xxt, _ := x.TimesDense(x.Transpose())
	this.XXt.Subtract(xxt)
	yxt, _ := y.TimesDense(x.Transpose())
	this.YXt.Subtract(yxt)	
}

func (this *KnownVarianceLRPosterior) GetSampler() func () *mx.DenseMatrix {
	if this.sampler == nil {
		PhiInv, err := this.Phi.Inverse()
		if err != nil { panic(err) }

		XXtpPhiInv, err := this.XXt.PlusDense(PhiInv)
		if err != nil { panic(err) }

		Omega, err := XXtpPhiInv.Inverse()
		if err != nil { panic(err) }

		MPhiInv, err := this.M.TimesDense(PhiInv)
		if err != nil { panic(err) }

		YXtpMPhiInv, err := this.YXt.PlusDense(MPhiInv)
		if err != nil { panic(err) }

		Mxy, err := YXtpMPhiInv.TimesDense(Omega)
		if err != nil { panic(err) }

		this.sampler = stat.MatrixNormal(Mxy, this.Sigma, Omega)
	}
	return this.sampler
}

func (this *KnownVarianceLRPosterior) Sample() (A *mx.DenseMatrix) {
	return this.GetSampler()()
}



/*
	If Y ~ N(AX, Sigma, I)
	and A ~ N(M, Sigma, Phi)
	this returns a sampler for P(A|X,Y,Sigma,M,Phi)
*/
func KnownVariancePosterior(Y, X, Sigma, M, Phi *mx.DenseMatrix) func() (A *mx.DenseMatrix) {
	o := Y.Rows()
	i := X.Rows()
	n := Y.Cols()
	if n != X.Cols() {
		panic("X and Y don't have the same number of columns")
	}
	if o != M.Rows() {
		panic("Y.Rows != M.Rows")
	}
	if i != M.Cols() {
		panic("Y.Rows != M.Cols")
	}
	if o != Sigma.Rows() {
		panic("Y.Rows != Sigma.Rows")
	}
	if Sigma.Cols() != Sigma.Rows() {
		panic("Sigma is not square")
	}
	if i != Phi.Rows() {
		panic("X.Rows != Phi.Rows")
	}
	if Phi.Cols() != Phi.Rows() {
		panic("Phi is not square")
	}

	Xt := X.Transpose()

	PhiInv, err := Phi.Inverse()
	if err != nil { panic(err) }

	XXt, err := X.TimesDense(Xt)
	if err != nil { panic(err) }

	XXtpPhiInv, err := XXt.PlusDense(PhiInv)
	if err != nil { panic(err) }

	Omega, err := XXtpPhiInv.Inverse()
	if err != nil { panic(err) }

	YXtpMPhiInv, err := Y.TimesDense(Xt)
	if err != nil { panic(err) }

	MPhiInv, err := M.TimesDense(PhiInv)
	if err != nil { panic(err) }

	err = YXtpMPhiInv.AddDense(MPhiInv)
	if err != nil { panic(err) }

	Mxy, err := YXtpMPhiInv.TimesDense(Omega)
	if err != nil { panic(err) }

	if false {
		fmt.Printf("Mxy:\n%v\n", Mxy)
		fmt.Printf("Sigma:\n%v\n", Sigma)
		fmt.Printf("Omega:\n%v\n", Omega)
	}

	return stat.MatrixNormal(Mxy, Sigma, Omega)
}
