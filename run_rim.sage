attach("ghost.py")
attach("halo_slopes.sage")

def test_on_rim(p,comp):
	g=ghost(p,1)
	v2 = slopes_on_rim(p))
	v1 = g.wadic_slopes(comp,num=len(v2))
	return v1 == v2

def write_to_file(filename,line):
	"""
	Writes the string line to the file ``filename``.

	Input:

	-	``filename`` -- string
	-	``line`` -- string
	"""
	f = open(filename,'a')
	f.write(line)
	f.close()

def run():
	for p in primes(3,100):
		for a in range(0,p-1,2):
			ans = test_on_rim(p,a)
			write_to_file("rim.data2",str(p)+","+str(a)+":"+str(ans)+"\n")
			print p,a,":",ans
			

