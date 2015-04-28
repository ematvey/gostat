package stat

import "testing"

func TestBetaInvCDF(t *testing.T) {
	dist := BetaInv_CDF(1, 4001)
	if dist(0.5) != 0.0001732284784849848 {
		t.Fail()
	}
}
