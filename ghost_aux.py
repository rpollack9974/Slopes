def separating_valuation(g,comp,j,z):
	"""returns the maximum valuation of z-z' where z' runs over all zeroes of the j-th coefficient of g on 
		component comp distinct from z"""
	zs = g.zeroes(comp,j)	
	return max([(z-ZZ(zp)).valuation(g.p) for zp in zs if zp!=z])+1

def total_separating_valuation(g,comp,j):
	"""returns v such that ord_p(z-z')<=v for all z,z' distinct zereos of j-th coefficient of g 
		on component comp"""
	return max([separating_valuation(g,comp,j,ZZ(z)) for z in g.zeroes(comp,j)])

def optimized_bound(g,comp,j,z):
	"""assuming all true zeroes of this coefficient are closer to a corresponding ghost zero (as expected)
		and that the ghost conjecture in predicts slopes correctly in lower indices, 
		returns a lower bound on the sum of valuations of the true zeroes near z to ensure that the j-th
		Newton point is high enough to lie at or above the semistable line UGH"""
	if j==1 and z==10:
		return 8
	else:
		dz = g.first_occurrence(comp,z)-1 ##classical dimension
		a = sum(g.slopes(z,num=dz,force_comp=comp))+(j-dz)*(z-2)/2
		b = 0
		p = g.p
		for zp in g.zeroes(comp,j):
			if zp != z:
				b += g.mult(comp,j,zp) * ((z-ZZ(zp)).valuation(p)+1)
		return (a-b)

def good_weight(g,comp,j,k):
	"""returns whether or not k is too close to a ghost zero of g in index j (too close = false)"""
	p = g.p
	good_weight = true
	for z in g.zeroes(comp,j):
		bnd = optimized_bound(g,comp,j,z)
		good_weight = ((ZZ(k)-z).valuation(p)+1 < bnd) and good_weight
	return good_weight


def available_test_weights(g,comp,i,verbose=false):
	"""assuming control over all lower indices j<i, returns a list of weights at which
		we should know exactly the valuation of a_i"""
	p = g.p
	vodd = [k for k in range(3*i,3*i+6) if k%2==1]
	veven = [k for k in range(3*i,3*i+6) if k%2==0]

	for k in range(3*i+6,6*i+3):
		if k % 2 == 1:
			#print "Trying",k
			dp = floor(k/3) - 1
			j = dp - i
			if good_weight(g,comp,j,k):
				vodd += [k]
			else:
				if verbose:
					print "eliminating",k,": problem at index",j

	for k in range(3*i+6,4*i+6):
		if k % 2 == 0:
			#print "Trying",k
			dp = floor(k/3) - 1
			j = dp - i
			if good_weight(g,comp,j,k):
				s = g.slopes(comp,k,force_comp=comp)
				if s[j-1] != s[j]:
					vodd += [k]
				else:
					if verbose:
						print "eliminating",k,": problem at index",j,"because of repeated slopes"
			else:
				if verbose:
					print "eliminating",k,": problem at index",j

	# vodd = [k for k in range(3*i,6*i+3) if k%2==1]
	# veven = [k for k in range(3*i,4*i+6) if k%2==0]

	## theta operator range
	vtheta = []
	for k in range(3,3*i):
		dk = floor(k/3)-1 
		if k == 3:
			dk = 0
		bool = true
		for j in range(1,dk+1):
			bool = bool and good_weight(g,0,j,k)
		kp = 2 - k
		for j in range(i-dk,i):
			bool = bool and good_weight(g,0,j,kp)
		if bool:
			v_cl = g.slopes(k,force_comp=0)[0:dk+1]  ## extra + 1 here for theta(eisen)
			v_kp = g.slopes(kp,force_comp=0)[0:i-1]
			y_cl = [sum(v_cl[0:a]) for a in range(0,len(v_cl)+1)]
			y_kp = [sum(v_kp[0:a]) for a in range(0,len(v_kp)+1)]
			y_theta_kp = [y_kp[a]+a*(k-1) for a in range(len(y_kp))]
			w = [y_cl[a] + y_theta_kp[i-a] for a in range(1,dk+2)]
			y_min = min(w)
			if verbose and w.count(y_min) != 1:
				print "multiple root in theta operator range for k=",k
			else:
				if y_min >= g.lambdainv(0,i) + (k-1)*i:
					if verbose:
						print "lower bound doesn't suffice in theta range for k=",k
				else:
					vtheta += [k]
					if verbose:
						print "For theta range, found for k=",k,"at",i,"-th slope ",y_min-sum(g.slopes(k,force_comp=0)[0:i-1])
		elif verbose:
			print k,"not a good weight"

	ans = vodd + veven

	ans = ans + vtheta
	ans.sort()

	return ans

def find_test_weights(g,comp,j,z):
	p = g.p
	v = separating_valuation(g,comp,j,z)
	tests = available_test_weights(g,comp,j)
	ks = [k for k in tests if (ZZ(k)-z).valuation(p) + 1 >= v]
	return ks

## implicitly assuming p=3, N=1
def locate_nearby_zero(g,comp,j,z):
	p = g.p
	ks = find_test_weights(g,comp,j,z)
	v = [(ZZ(z)-ZZ(k)).valuation(p) for k in ks]
	m = max(v)
	t = v.index(m)
	k = ks[t]

	dk = dimension_cusp_forms(1,k)
	d = floor(k/3)-1 ## this equals the dimenion of S_k(Gamma_1(3),omega^k)
	## When j matches exactly the dimension or the dimension plus 1 we know exactly
	## the height of the j-th Newton point (which is a break)
	if (j == d) or (j == d+1):
		## when k is odd, the symmetry of the slopes gives the height if j==d and we must
		## add k-1 more when j==d+1 for theta of Eisenstein series.  
		if k%2 == 1:
			y = (k-1)*d/2
		else:
			y = (k-1)*dk+(dimension_cusp_forms(g.p,k)-2*dk)*(k-2)/2
		if j == d+1:
			y = y + k-1
	elif j < d:
		jp = d-j
		y = sum([((ZZ(zz)-ZZ(k)).valuation(p)+1)*g.mult(comp,jp,zz) for zz in g.zeroes(comp,jp)])
		if k%2 == 0:
			y = y + (dk-jp)*(k-1) + (dimension_cusp_forms(g.p,k)-2*dk) * (k-2)/2
		else:
			y = y + (j-d/2)*(k-1)
	else:
		y = sum(g.slopes(k,force_comp=0)[0:j])
	dz = g.first_occurrence(comp,z)-1 ##classical dimension
	y2 = sum(g.slopes(z,num=dz,force_comp=comp))+(j-dz)*(z-2)/2
	if y2 > y:
		print "For ghost zero",z,"using k =",k,"found zero closer than v=",(ZZ(k)-ZZ(z)).valuation(p)+1
		print "  			v(a_j(",k,"))=",y," and v(a_j(",z,"))>=",y2
		return k
	else:
		print "Weight",k,"gave too high a valuation"


#****JUNKY CODE BELOW

def locate_nearby_zeroes(g,comp,j):
	for z in g.zeroes(comp,j):
		locate_nearby_zero(g,comp,j,z)

def available_test_weights_optimistic(g,comp,i):
	v = range(3,6*i+2,2) + range(4,4*i+5,2)
	v.sort()
	return v

def nearby_test_weight(g,comp,j,z):
	p = g.p
	ks = available_test_weights_optimistic(g,comp,j)
	v = separating_valuation(g,comp,j,z)
	w = [(ZZ(z)-ZZ(k)).valuation(p) for k in ks]
	m = max(w)
	t = w.index(m)
	k = ks[t]
	assert (ZZ(z)-ZZ(k)).valuation(p) + 1 >= v,"oops"
	return k

def locate_nearby_zero_optimistic(g,comp,j,z):
	t = nearby_test_weight(g,comp,j,z)
	y = sum(g.slopes(t,num=j,force_comp=0))
	dz = dimension_cusp_forms(1,z)
	y2 = sum(g.slopes(z,num=dz,force_comp=comp))+(j-dz)*(z-2)/2
	if y2 > y:
		return true,t
	else:
		return false,t








