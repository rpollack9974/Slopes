def do_it(N):
	v = []
	for j in range(0,N):
		ans =0
		for z in g.zeroes(0,j):
			if z%5!=0:
				ans += g.mult(0,j,z)
			else:
				ans += g.mult(0,j,z)*1.5
		v += [ans]
	return v

def bn(n):
	R=PolynomialRing(QQ,'x')
	if n==0:
		return 1/2+0*x
	elif n==1:
		return 1 + 0*x
	else:
		return 1/21*(-(x-7)*(64*x-7)*bn(n-1).derivative() + (32*(n-1)*x-56*(n-1)+42)*bn(n-1)-2*(n-1)*(2*(n-1)-1)*(11*x+7)*bn(n-2))
