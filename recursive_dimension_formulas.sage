class rbdata(SageObject):
	def __init__(self,p,krbar,split,mult):
		self.p = p
		self.krbar = krbar
		self.split = split
		self.mult = mult

	def __repr__(self):
		s1 = "Mod "+str(self.p)+" repn with Serre weight "+str(self.krbar)+", "
		if self.split:
			s2 = "split at p, "
		else:
			s2 = "non-split at p, "

		return repr(s1 + s2 + "mults: "+str(self.mult))

## return unique integer r(k) with r(k) = k (mod N) and 2 < r(k) <= N+2
def r(k,N):
	return (k-3)%(N)+3

##Returns the change of the rhobar-dimension when k is increased by p^2-1
##(see Corollary 5.9 in currect draft)
def Q_dt(rb):
	p = rb.p
	if not rb.split:
		if rb.krbar%(p-1)!=2%(p-1):
			return 2*rb.mult[0]
		else:
			return rb.mult[0]+rb.mult[1]
	else:
		if rb.krbar%(p-1)!=2%(p-1):
			return 2*rb.mult[0]+2*rb.mult[2]
		else:
			return rb.mult[0]+rb.mult[1]+2*rb.mult[2]

def S(k,t,rbdata):
	p = rbdata.p
	t = t%(p-1)
	if k % (p-1) != (rbdata.krbar + 2*t) % (p-1):
		return 0
	elif k>=p^2+1:
		if k % (p^2-1) >= 2:
			m = floor(k/(p^2-1))
		else:
			m = floor(k/(p^2-1))-1
		return Q_dt(rbdata)*m+S(k-m*(p^2-1),t,rbdata)
	elif k >= p+3:
		return S(k-(p+1),t-1,rbdata) + S(r(k,p-1),t,rbdata) + S(p+3-r(k,p-1),t-k+2,rbdata)
	elif k == p+2:
		return S(3,t,rbdata) + S(p,t-1,rbdata)
	elif rbdata.split == False:
		if rbdata.krbar % (p-1) != 2 % (p-1):
			if k == rbdata.krbar and t==0:
				return rbdata.mult[0]
			else:
				return 0
		else:
			if t == 0 and k == 2:
				return rbdata.mult[0]
			elif t == 0 and k == p+1:
				return rbdata.mult[1]
			else:
				return 0
	else:
		if rbdata.krbar % (p-1) != 2 % (p-1):
			if k == rbdata.krbar and t == 0:
				return rbdata.mult[0]
			elif k == p+1-rbdata.krbar and t == p-rbdata.krbar:
				return rbdata.mult[2]
			else:
				return 0
		else:
			if k == rbdata.krbar and t == 0:
				return rbdata.mult[0]
			elif k == p+1-rbdata.krbar and t == p-rbdata.krbar:
				return rbdata.mult[2]
			elif k == p+1 and t == 0:
				return rbdata.mult[1]
			else:
				return 0

# def Sp(k,t,rbdata):
# 	p = rbdata.p
# 	t = t%(p-1)
# 	if k % (p-1) != (rbdata.krbar + 2*t) % (p-1):
# 		return 0
# 	if k<p+1:
# 		kzt = (rbdata.krbar + 2*t - 2)%(p-1)+2
# 		if rbdata.split == False and rbdata.krbar % (p-1) != 2:
# 			if t <= kzt-2:
# 				return 2*rbdata.mult[0]
# 			else:
# 				return 0
# 		else:
# 			tot = 0
# 			for j in range(0,kzt-1):
# 				if (t+j-(kzt-2)) % (p-1) == 0:
# 					if rbdata.krbar % (p-1) == 2:
# 						tot += rbdata.mult[2]
# 					else:
# 						tot += rbdata.mult[0]
# 				if (t-j) % (p-1) == 0:
# 					tot += rbdata.mult[0]
# 				if rbdata.split:
# 					if (t+j-(kzt-2)) % (p-1) == (p - rbdata.krbar) % (p-1):
# 						tot += 1
# 					if (t-j) % (p-1) == (p - rbdata.krbar) % (p-1):
# 						tot += 1
# 			return tot
# 	elif not rbdata.split:
# 		s = num_steps(k,p)
# 		if rbdata.krbar % (p-1) != 2:
# 			#return Sp(k-(p-1),t,rbdata) + 2*rbdata.mult[0]
# 			return 2*rbdata.mult[0]*s + Sp((k-2)%(p-1)+2,t,rbdata)
# 		else:
# 			#return Sp(k-(p-1),t,rbdata) + rbdata.mult[0] + rbdata.mult[1]
# 			return (rbdata.mult[0] + rbdata.mult[1])*s + Sp((k-2)%(p-1)+2,t,rbdata)
# 	else:
# 		if rbdata.krbar % (p-1) != 2:
# 			return Sp(k-(p-1),t,rbdata) + 2*rbdata.mult[0] + 2*rbdata.mult[2]
# 		else:
# 			return Sp(k-(p-1),t,rbdata) + rbdata.mult[0] + rbdata.mult[1] + 2*rbdata.mult[2]

def Sp(k,t,rbdata):
	p = rbdata.p
	t = t%(p-1)
	if k % (p-1) != (rbdata.krbar + 2*t) % (p-1):
		return 0
	elif k>=p+1:
		if k % (p-1) >= 2:
			m = floor(k/(p-1))
		else:
			m = floor(k/(p-1))-1
		return Q_dt(rbdata)*m+Sp(k-m*(p-1),t,rbdata)
	else:
		kzt = (rbdata.krbar + 2*t - 2)%(p-1)+2
		if rbdata.split == False and rbdata.krbar % (p-1) != 2 % (p-1):
			if t <= kzt-2:
				return 2*rbdata.mult[0]
			else:
				return 0
		else:
			tot = 0
			for j in range(0,kzt-1):
				if (t+j-(kzt-2)) % (p-1) == 0:
					if rbdata.krbar % (p-1) == 2 % (p-1):
						tot += rbdata.mult[1]
					else:
						tot += rbdata.mult[0]
				if (t-j) % (p-1) == 0:
					tot += rbdata.mult[0]
				if rbdata.split:
					if (t+j-(kzt-2)) % (p-1) == (p - rbdata.krbar) % (p-1):
						tot += rbdata.mult[2]
					if (t-j) % (p-1) == (p - rbdata.krbar) % (p-1):
						tot += rbdata.mult[2]
			return tot

def num_steps(k,p):
	r = k % (p-1)
	q = (k-r)/(p-1)
	if r != 0:
		return q
	else:
		return q-1
