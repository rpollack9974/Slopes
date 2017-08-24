def test_distinct_wadic_slopes_fixed_rhobar(rb,verbose=false):
	p = rb.p
	M = Infinity
	for t in range(0,p-1):
		g = ghost(rbdata=rb,twist=t,terms=4*rb.mult[0]*p)
		v = g.wadic_slopes()
		w = v[2*rb.mult[0]*p:4*rb.mult[0]*p]
		y = [w[a]-v[a] for a in range(len(w))]
		if min(y)!=max(y):
			print "AP failed!!:",t
		w = [v[a]-v[a-1] for a in range(1,len(v))]
		if verbose:
			print t,min(w)
		M = min(M,min(w))
	return M

def test_distinct_wadic_slopes(p,split):
	M = Infinity
	for krbar in range(4,p+2,2):
		rb = rbdata(p,krbar,split,[1,1,1])
		print "krbar =",krbar,": min difference in w-adic slopes is:",test_distinct_wadic_slopes_fixed_rhobar(rb)
	return "Done"
