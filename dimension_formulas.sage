## Returns the dimension of S_k(Gamma_0(N))_{rhobar otimes omega^t} where
## (locally reducible) rhobar occurs in S_{kt} with kt < p+1 with multiplicity m.
## THESE FORMULAS SHOULD BE CHECKED!!!  Almost certainly don't work for kt=2
def levelN_rhobar_dimension_reducible(rbdata,k,t=0):
	p=rbdata[0]
	k0=rbdata[1]
	mult=rbdata[2]
	if k%(p-1)!=(k0+2*t)%(p-1):
		return 0
	else:
		if t<p+1-k0:
			n=(k-k0-t*(p+1))/(p-1)
			return (2*floor(n/(p+1))+d1(rbdata,n))*mult
		else:
			tp=t-(p+1-k0)
			k0p=2*p-k0+2
			n=(k-k0p-tp*(p+1))/(p-1)
			return (2*floor(n/(p+1))+d2(rbdata,n))*mult

	# p=rbdata[0]
	# kt=rbdata[1]
	# m=rbdata[2]
	# t=t%(p-1)
	# if k%(p-1)!=(kt+(p+1)*t)%(p-1):
	# 	return 0

	# if t<p-kt+1:
	# 	k=k-t*(p+1)
	# 	j=(k-2)%(p-1)+2  # so 2 <= j < p+1
	# 	m=(k-j)/(p-1)
	# 	mm=m%(p+1)
	# 	if mm<kt:
	# 		return 1+2*floor(m/(p+1))
	# 	else:
	# 		return 2+2*floor(m/(p+1))
	# else:
	# 	t=t-(p-kt+1)
	# 	k=k-t*(p+1)
	# 	j=2*p-kt+2. ## This is the weight after the drop in the theta cycle
	# 	#j=(k-2)%(p-1)+2  # so 2 <= j < p+1
	# 	m=(k-j)/(p-1)
	# 	mm=m%(p+1)
	# 	if (mm<p-kt+1) and (mm>=0):
	# 		return 1+2*floor(m/(p+1))
	# 	else:
	# 		return 2+2*floor(m/(p+1))


def levelNp_rhobar_dimension_reducible(rbdata,k,t=0):
	p=rbdata[0]
	k0=rbdata[1]
	mult=rbdata[2]
	if k%(p-1)!=(k0+(p+1)*t)%(p-1):
		return 0
	k0t=(k-2)%(p-1)+2  # so 2 <= k0t < p+1
	m = (k-k0t)/(p-1)
	if k0t < t+2:  
		m=m-1
	return 2*(m+1)*mult


def levelNp_pnew_rhobar_dimension_reducible(rbdata,k,t=0):
	return levelNp_rhobar_dimension_reducible(rbdata,k,t)-2*levelN_rhobar_dimension_reducible(rbdata,k,t)




def r(n,m):
	"""Returns n mod m normalized so that 0<=r(n,m)<m"""
	return n%m

def d1(rbdata,n):
	"""Returns 1 if n<k0 and 2 otherwise"""
	p=rbdata[0]
	k0=rbdata[1]
	if r(n,p+1)<k0:
		return 1
	else:
		return 2

def d2(rbdata,n):
	"""Returns 1 if n<p+1-k0 and 2 otherwise"""
	p=rbdata[0]
	k0=rbdata[1]
	if r(n,p+1)<p+1-k0:
		return 1
	else:
		return 2


def lzp(rbdata,i):
	p=rbdata[0]
	k0=rbdata[1]
	return k0+(i-1)*(p-1)

def hzm(rbdata,i):
	p=rbdata[0]
	k0=rbdata[1]
	return k0+(i-2)*(p-1)

def hzp(rbdata,i):
	p=rbdata[0]
	k0=rbdata[1]
	return k0+floor((i-2)/2)*(p^2-1)+a(rbdata,i)

def a(rbdata,i):
	p=rbdata[0]
	k0=rbdata[1]
	if i%2==0:
		return (p-1)*(k0-1)
	else:
		return (p-1)*(p)

def lzm(rbdata,i):
	p=rbdata[0]
	k0=rbdata[1]
	return k0+floor((i-1)/(2*p))*(p^2-1)+b(rbdata,i)

def b(rbdata,i):
	p=rbdata[0]
	k0=rbdata[1]
	i=(i-1)%(2*p)+1 #So 1<=i<=2p
	if i<=2*k0:
		return floor(i/2)*(p-1)
	else:
		return floor((i+1)/2)*(p-1)

def ws(rbdata,i):
	p=rbdata[0]
	k0=rbdata[1]
	return (p+1)*(floor((i-2)/2)+floor((i-1)/(2*p)))-2*i+3+(a(rbdata,i)+b(rbdata,i))/(p-1)

def f(rbdata,i,t):
	p=rbdata[0]
	k0=rbdata[1]
	i=(i-1)%(2*p)+1
	if i%2==0:
		if i<=2*t:
			return p-1
		else:
			return 0
	else:
		if i-p<2*t and i-p>=0:	
			return p-1
		else:
			return 0

def lzmt(rbdata,i,t):
	p=rbdata[0]
	k0=rbdata[1]
	if t<p+1-k0:

		return k0+2*t+floor((i-1)/(2*p))*(p^2-1)+b(rbdata,i)-f(rbdata,i,t)
	else:
		return k0+2*t+floor((i-1)/(2*p))*(p^2-1)+b(rbdata,i)-f(rbdata,i,t-p+k0)-(p-1)

def k0t(rbdata,t):
	p=rbdata[0]
	k0=rbdata[1]
	return ((k0+2*t)-2)%(p-1)+2

def lzmt2(rbdata,i,t):
	p=rbdata[0]
	k0=rbdata[1]
	r=(i)%(2*p)
	C=ceil((r-1)/2)
	if t<(2*k0-1):
		if i%2==0:
			if C>t:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*C
			else:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C-1)
		else:
			if ((C-t)%(p-1)<k0):
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*C
			else:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
	else:
		if i%2==0:
			if C>t:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
			else:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C)
		else:
			if ((C-t)%(p-1)<k0):
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
			else:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+2)

def lzmt3(rbdata,i,t):
	p=rbdata[0]
	k0=rbdata[1]
	r=(i)%(2*p)
	C=ceil((r-1)/2)
	if i%2==1:
		if k0t(rbdata,t)>=t+2:
			if C>=t and C-t<k0 and t<p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C)
			elif C>=t and C-t<k0 and t>=p+1-k0:
	 			return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C)
			elif C>=t and C-t>=k0 and t<p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
			elif C<t and p+1+C-t>=k0 and t<p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C)
			elif C<t and p+1+C-t<k0 and t<p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C)
			elif C<t and p+1+C-t<k0 and t>=p+1-k0:
				if r<=2*k0+1:
					return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C)
				else:
					return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
			elif C>=t and C-t>=k0 and t>=p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C)
			elif C<t and p+1+C-t>=k0 and t>=p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
		else:
			if C>=t and C-t<k0 and t<p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
			elif C>=t and C-t<k0 and t>=p+1-k0:
	 			return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+2)
			elif C>=t and C-t>=k0 and t<p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+2)
			elif C<t and p+1+C-t>=k0 and t<p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
			elif C<t and p+1+C-t<k0 and t<p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
			elif C<t and p+1+C-t<k0 and t>=p+1-k0:
				if r<=2*k0+1:
					return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
				else:
					return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+2)
			elif C>=t and C-t>=k0 and t>=p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)
			elif C<t and p+1+C-t>=k0 and t>=p+1-k0:
				return k0t(rbdata,t)+floor((i)/(2*p))*(p^2-1)+(p-1)*(C+1)+1/2
	else:
		return -1000 



