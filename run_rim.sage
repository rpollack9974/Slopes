attach("ghost.py")
attach("lambda.sage")

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
			write_to_file("rim.data",str(p)+","+str(a)+":"+str(ans)+"\n")
			print p,a,":",ans
			

