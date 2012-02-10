package stat

import (
	"fmt"
	mtx "gomatrix.googlecode.com/hg/matrix"
	"testing"
)

// Covariance test against R
func TestCovariance(t *testing.T) {
	fmt.Println("Covariance test against R")
	data := GetData()
	out := SCov(data)

	//known values
	dist := [...]float64{1.773643,-0.3381504,0.4803343,-0.8050336,0.3154475,0.2026154,0.5576372,-0.3982494,0.2083944,0.1306447,
-0.3381504,0.8108673,0.1696915,0.5268912,-0.1279332,-0.0810701,0.4332097,-0.01344127,0.3683854,-0.1976432,
0.4803343,0.1696915,0.4986279,0.02131043,0.2643708,-0.3654575,0.6973307,-0.3021538,-0.221391,-0.07286902,
-0.8050336,0.5268912,0.02131043,0.7320278,0.01051328,-0.4524992,0.0321159,-0.1932133,0.1346962,0.003491767,
0.3154475,-0.1279332,0.2643708,0.01051328,0.4088249,-0.3923959,0.2455576,-0.3880454,-0.2810396,0.1249063,
0.2026154,-0.0810701,-0.3654575,-0.4524992,-0.3923959,0.8562731,-0.3241916,0.5470234,0.3595004,-0.1266668,
0.5576372,0.4332097,0.6973307,0.0321159,0.2455576,-0.3241916,1.149185,-0.1579107,-0.3901415,-0.3014224,
-0.3982494,-0.01344127,-0.3021538,-0.1932133,-0.3880454,0.5470234,-0.1579107,0.7091816,-0.2188541,-0.3055508,
0.2083944,0.3683854,-0.221391,0.1346962,-0.2810396,0.3595004,-0.3901415,-0.2188541,1.283379,0.2383891,
0.1306447,-0.1976432,-0.07286902,0.003491767,0.1249063,-0.1266668,-0.3014224,-0.3055508,0.2383891,0.2644874}

	cols := data.Cols()
	known := mtx.Zeros(cols, cols)
	for i := 0; i < cols; i++ {
		for j := i + 1; j < cols; j++ {
			known.Set(i, j, dist[i*cols+j])
		}
	}

	// check
	for i := 0; i < cols; i++ {
		for j := i + 1; j < cols; j++ {
			if !check(out.Get(i, j), known.Get(i, j)) {
				t.Error()
				fmt.Println(i, j, out.Get(i, j), known.Get(i, j))
			}
		}
	}
}

// Variance test against R
func TestVariance(t *testing.T) {
	fmt.Println("Variance test against R")
	data := GetData()
	out := SVar(data)

	//known values
	dist := [...]float64{1.7736433,0.8108673,0.4986279,0.7320278,0.4088249,0.8562731,1.1491851,0.7091816,1.2833793,0.2644874}





	cols := data.Cols()
				fmt.Println("cols: ", cols)

	known := NewVector(cols)
	for i := 0; i < cols; i++ {
			known.Set(i, dist[i])
	}

	// check
	for i := 0; i < cols; i++ {
			if !check(out.Get(i), known.Get(i)) {
				t.Error()
			}
				fmt.Println(i, out.Get(i), known.Get(i))

	}
}


