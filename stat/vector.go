package stat

type Vector struct {
	A []float64 // data
	L int       // length
}

func NewVector(length int) (v *Vector) {
	v = new(Vector)
	v.L = length
	v.A = make([]float64, length)
	return v
}

func (v Vector) Set(i int, x float64) {
	v.A[i] = x
}

func (v Vector) Get(i int) float64 {
	return v.A[i]
}

func (v Vector) Swap(i int, j int) {
	x := v.A[i]
	v.A[i] = v.A[j]
	v.A[j] = x
}

func (v Vector) Len() int {
	return v.L
}


