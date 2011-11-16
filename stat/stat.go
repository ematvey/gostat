//Extended set of distributions.
package stat

import (
	"math/rand"
	"time"
)

var Seed func(int64) = rand.Seed

func TimeSeed() {
	Seed(time.Nanoseconds())
}
