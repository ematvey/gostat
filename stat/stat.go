//Extended set of distributions.
package stat

import (
	"rand"
	"time"
)

var Seed func(int64) = rand.Seed

func TimeSeed() {
	Seed(time.Nanoseconds())
}
