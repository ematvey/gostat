// test for F_Qtl_For
package stat

import (
	"testing"
	"fmt"
)

// test against known values
func TestF_Qtl_For(t *testing.T) {
	fmt.Println("test for F_Qtl_For")
	var df1, df2, x, y, p float64
	df1=3
	df2=3
	x=0.46
	cdf:=F_CDF(df1, df2)
	p=cdf(x)
	y = F_Qtl_For(df1, df2, p)

	if !check(x, y){
		t.Error()
	}
}


