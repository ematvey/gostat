package bayes

import (
	"testing"
	"fmt"
)

// test fmin against known values

func TestMin(t *testing.T) {
	x := fmin(parab, 0, 1, 1e-6)
	y := 0.3333333333333333333
			if !check(x, y){
				t.Error()
				fmt.Println(x, y)

			}
}

