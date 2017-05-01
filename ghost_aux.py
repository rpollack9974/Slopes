def separating_valuation(g,j,comp,z):
	"""returns the maximum valuation of z-z' where z' runs over all zeroes of the j-th coefficient of g on 
		component comp distinct from z"""
	zs = g.zeroes(j,comp)	
	return max([(z-zp).valuation(g.p) for zp in zs if zp!=z])+1

def total_separating_valuation(g,j,comp):
	"""returns v such that ord_p(z-z')<=v for all z,z' distinct zereos of j-th coefficient of g 
		on component comp"""
	return max([separating_valuation(g,j,comp,ZZ(z)) for z in g.zeroes(j,comp)])

def optimized_bound(g,j,comp,z):
	"""assuming all true zeroes of this coefficient are closer to a corresponding ghost zero (as expected)
		and that the ghost conjecture in predicts slopes correctly in lower indices, 
		returns a lower bound on the sum of valuations of the true zeroes near z to ensure that the j-th
		Newton point is high enough to lie at or above the semistable line UGH"""
	dz = g.first_occurrence(comp,z)-1 ##classical dimension
	a = sum(g.slopes(z,num=dz,force_comp=comp))+(j-dz)*(z-2)/2
	b = 0
	p = g.p
	for zp in g.zeroes(j,comp):
		if zp != z:
			b += g.mult(j,comp,zp) * ((z-zp).valuation(p)+1)
	return (a-b)


