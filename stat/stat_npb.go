package stat

func CRP_PDF(α float64)  func(x []int64) float64 {
	return func(x []int64) float64 {
		n := int64(len(x));
		counts := make([]int64, int(α*log(float64(len(x)))));
		sum := fZero;
		
		p := fOne;
		
		for i:=iZero; i<n; i++ {
			if x[i] >= int64(len(counts)) {
				counts = copyInt64(counts, 2*int64(len(counts)))
			}
			
			if counts[x[i]] == 0 {
				p *= α/(sum+α);
			} else {
				p *= float64(x[i])/(sum+α);
			}
			
			counts[x[i]] += 1;
			sum += 1;
		}
		return p;
	}
}
func CRP_LnPDF(α float64) func(x []int64) float64 {
	return func(x []int64) float64 {
		n := int64(len(x));
		counts := make([]int64, int(α*log(float64(len(x)))));
		sum := fZero;
		
		p := fZero;
		
		for i:=iZero; i<n; i++ {
			if x[i] >= int64(len(counts)) {
				counts = copyInt64(counts, 2*int64(len(counts)))
			}
			
			if counts[x[i]] == 0 {
				p += log(α)-log(sum+α);
			} else {
				p += log(float64(x[i]))-log(sum+α);
			}
			
			counts[x[i]] += 1;
			sum += 1;
		}
		return p;
	}
}
