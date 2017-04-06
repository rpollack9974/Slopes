attach("ghost.py")
attach("lambda.sage")

for p in primes(3,100):
	for a in range(0,p-1,2):
		print p,a,":",test_on_rim(p,1,a)

