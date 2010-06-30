package stat

import (
	"rand";
	"time"
)

var Seed func(int64) = rand.Seed

func init() {
	Seed(time.Nanoseconds());
}
