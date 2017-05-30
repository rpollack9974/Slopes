class ghost_proof(ghost):
	def __init__(self,p,N,num_coefs=10):
		"""Initializes the ghost series

		Inputs:
			p - prime
			N - tame level
			num_coefs (optional) - number of coefficients computed
		"""
		self.p=p
		self.N=N
		self.compute_ghost_series(num_coefs=num_coefs)
		self.zero_data = {'proven_info':True}
		for comp in range(self.p-1):
			for j in range(num_coefs+1):
				for z in self.zeroes(comp,j):
					## Putting 1 down as a lower bound for the distance between each ghost
					## zero and a corresponding true zero is simply the statement that all
					## of the zeroes of the a_i are in the central disc (i.e. Roe's result)
					self.zero_data[(comp,j,z)] = [[1 for a in range(self.mult(comp,j,z))],self.mult(comp,j,z)]

	def sep_ghost_zero(self,comp,j,z):
		"""returns valuation v so that all other ghost zeroes z' on component comp and index j 
		satisfy v(w_z-w_z')<=v"""
		p = self.p
		assert self.zeroes(comp,j).count(z)>0,"Not a ghost zero in sep_ghost_zero"
		zs = self.zeroes(comp,j)
		return max([(ZZ(zp)-ZZ(z)).valuation(p)+1 for zp in zs if zp != z])

	def sep_nth_level_ghost_zero(self,comp,j,z,n):
		"""returns valuation v so that all other ghost zeroes z' on component comp and index j 
		which have multiplicity >= n satisfy v(w_z-w_z')<=v"""
		p = self.p
		assert self.mult(0,j,z)>=n,"Not an n-th level ghost zero in sep_nth_level_ghost_zero"
		zs = self.nth_level_zeroes(comp,j,n)
		return max([(ZZ(zp)-ZZ(z)).valuation(p)+1 for zp in zs if zp != z])

	def sep_all_ghost_zeroes(self,comp,j):
		return max([self.sep_ghost_zero(comp,j,z) for z in self.zeroes(comp,j)])

	def close_to_ghost_zero(self,comp,j,k,z):
		p = self.p
		return (ZZ(k)-ZZ(z)).valuation(p)+1 >= min(self.zero_data[(comp,j,z)][0])

	def close_to_nth_level_ghost_zero(self,comp,j,k,z,n):
		p = self.p
		return (ZZ(k)-ZZ(z)).valuation(p)+1 >= self.zero_data[(comp,j,z)][0][n-1]

	def good_weight(self,comp,j,k,verbose=false):
		"""returns whether or not the bounds recorded in self.zero_data are sufficient to bound
		k away from all ghost zeroes on component comp and index j"""
		p = self.p
		bool = true
		for z in self.zeroes(comp,j):
			if self.zero_data.has_key((comp,j,z)):
				bool2 = (ZZ(k)-ZZ(z)).valuation(p)+1 < min(self.zero_data[(comp,j,z)][0])
				bool = bool and bool2
				if not bool2 and verbose:
					print "Too close to",z
			else:
				return False
		return bool

	def y_lowerbnd(self,comp,j,k):
		p = self.p
		ans = 0
		for kss in self.zeroes(comp,j):
			if not self.close_to_ghost_zero(comp,j,k,kss):
				ans += ((ZZ(kss)-k).valuation(p)+1)*self.mult(comp,j,kss)
			else:
				for bnd in self.zero_data[(comp,j,kss)][0]:
					if (ZZ(kss)-k).valuation(p)+1 < bnd:
						ans += (ZZ(kss)-k).valuation(p)+1
					else:
						ans += bnd
		return ans

	def known_yvals(self,comp,k):
		v = []
		for j in range(self.num_coefs(comp)):
			ans_bool, ans_val = self.prove_yval(comp,j,k)
			if ans_bool:
				v += [ans_val]
			else:
				v += ["?"]
		return v

	def prove_yval(self,comp,j,k,verbose=False):
		"""tries to compute the exact valuation of P_j(w_k) on component comp"""
		### NEED TO GO BACK AND FIX FACT THAT COEFFICIENTS ARE ONLY ACCURATE MODULO NON_CLASSICAL PART
		p = self.p
		### if we already known the answer, it is returned
		if self.good_weight(comp,j,k):
			if verbose:
				print "We can determine the y-value at",k,"exactly because it is a good weight for a_",j
			return True,self.y(comp,j,k)
		elif k == 3*j or k == 3*j+1 or k == 3*j+2:
			##We are then at the index of the theta of Eisenstein
			if verbose:
				print "We can determine the y-value at",k,"exactly because it is the index of theta(Eisen)"
			if k%2 == 1:
				return True,(k-1)*(j-1)/2 + (k-1)
			else:
				d = dimension_cusp_forms(1,k)
				return True,d*(k-1) + (j-2*d-1)*(k-2)/2 + (k-1)
		elif k == 3*j+3 or k == 3*j+4 or k == 3*j+5:
			##We are then at the index of the last classical form and symmetry gives us the answer
			if verbose:
				print "We can determine the y-value at",k,"exactly because it is the index of last classical form"
			if k%2 == 1:
				return True,(k-1)*(j)/2
			else:
				d = dimension_cusp_forms(1,k)
				return True,d*(k-1) + (j-2*d)*(k-2)/2
		elif k%2 == 1 and k < 6*j+3 and k > 3*j:
			##We are then in an odd weight where symmetry can hopefully reflex us to a point
			##where we know the answer
			dp = floor(k/3)-1
			jj = dp - j
			if self.good_weight(comp,jj,k):
				if verbose:
					print "We can determine the y-value at",k,"exactly by using AL-symmetries"
				return True,self.y(comp,jj,k) + (dp/2-jj)*(k-1)
			else:
				return False,"not a good weight"
		elif k%2 == 0 and k < 4*j+6 and k > 3*j:
			##We are then in an even weight where symmetry can hopefully give us the answer
			dp = floor(k/3)-1
			d = dimension_cusp_forms(1,k)
			ptame = []
			## computes valuations coming from level N
			for i in range(d+1):
				if self.good_weight(comp,i,k):
					v = self.y(comp,i,k)
					ptame += [(v,"=")]
				else:
					v = self.y_lowerbnd(comp,i,k)
					ptame += [(v,">=")]
			## uses symmetry of p-refinement to compute remaining valuations
			for i in range(d+1,2*d+1):
				ptame += [(ptame[2*d-i][0]+(k-1)*(i-d),ptame[2*d-i][1])]
			dnew = dp - 2*d 
			## uses explicit formula for valuations of coefficients arising from p-new forms
			Pnew = Upnew_char_poly(k)
			Pnew_coefs = Pnew.padded_list()
			pnew = [(Pnew_coefs[a].valuation(p),"=") for a in range(Pnew.degree()+1)]
			#print ptame
			#print pnew
			## forms the product of forms from level N and p-new forms
			ans = []
			for a in range(len(ptame)):
				b = j - a
				if b < len(pnew):
					ans += [(ptame[a][0] + pnew[b][0],ptame[a][1])]
			#print ans
			ans_vals = [ans[a][0] for a in range(len(ans))]
			m = min(ans_vals)
			ind = ans_vals.index(m)
			if ans_vals.count(m)==1 and ans[ind][1] == "=":
				if verbose:
					print "We can determine the y-value at",k,"exactly by using oldform symmetries"
				return True,m
			else:
				if ans_vals.count(m)>1:
					return False,"min not unique"
				else:
					return False,"not a good weight"
		elif k < 3*j:
			d = floor(k/3) - 1 + 1 #+1 for theta(E_k)
			pcl = []
			## computes valuations coming from classical forms and theta of eisen
			for i in range(d+1):
				if self.good_weight(comp,i,k):
					v = self.y(comp,i,k)
					pcl += [(v,"=")]
				else:
					v = self.y_lowerbnd(comp,i,k)
					pcl += [(v,">=")]
			ptheta = []
			for i in range(j+1):
				if self.good_weight(comp,i,2-k):
					v = self.y(comp,i,2-k)
					ptheta += [(v,"=")]
				else:
					v = self.y_lowerbnd(comp,i,2-k)
					ptheta += [(v,">=")]
			print pcl
			print ptheta
			ptheta = [(ptheta[a][0] + a*(k-1),ptheta[a][1]) for a in range(len(ptheta))]
			print ptheta
			## forms the product of classical forms and forms in image of theta
			ans = []
			for a in range(len(pcl)):
				b = j - a
				v = pcl[a][0] + ptheta[b][0]
				if pcl[a][1] == ">=" or ptheta[b][1] == ">=":
					ans += [(v,">=")]
				else:
					ans += [(v,"=")]
			print ans
			ans_vals = [ans[a][0] for a in range(len(ans))]
			m = min(ans_vals)
			ind = ans_vals.index(m)
			if ans_vals.count(m)==1 and ans[ind][1] == "=":
				if verbose:
					print "We can determine the y-value at",k,"exactly via the theta operator"
				return True,m
			else:
				if ans_vals.count(m)>1:
					return False,"min not unique"
				else:
					return False,"theta fail"
		else:
			return False,"didn't work"

	def available_sample_weights(self,comp,j):
		return [k for k in range(3,6*j+6) if self.prove_yval(0,j,k)[0]]

	def find_first_level_true_zero(self,comp,j,z,verbose=False):
		assert self.zeroes(comp,j).count(z) > 0,"Not a ghost zero"
		if verbose:
			print "Attempting to find zero near",z
		r = self.sep_ghost_zero(comp,j,z)
		ks = self.available_sample_weights(comp,j)
		p = self.p
		v = [(ZZ(k)-ZZ(z)).valuation(p)+1 for k in ks]
		m = max(v)
		while m>=r:
			ind = v.index(m)
			k = ks[ind]
			if verbose:
				print "Trying sample weight",k
			sample = self.prove_yval(comp,j,k,verbose=verbose)
			if sample[0]:
				d = dimension_cusp_forms(1,z)
				bnd = self.y_lowerbnd(comp,d,z) + (z-2)/2 * (j-d)
				if verbose:
					print "Valuation at",k,"is",sample[1],"while valuation at",z,"is at least",bnd
				if sample[1] < bnd:
					if verbose:
						print "Success.  There is a true zero closer to",z,"than to",k,"\n"
					self.zero_data[(comp,j,z)][0][0] = (ZZ(k)-ZZ(z)).valuation(p) + 1.0001
					return True
				else:
					ks[ind] = z+1
					v = [(ZZ(k)-ZZ(z)).valuation(p)+1 for k in ks]
					m = max(v)					
			else:
				ks[ind] = z+1
				v = [(ZZ(k)-ZZ(z)).valuation(p)+1 for k in ks]
				m = max(v)

		if verbose:
			print "Failed\n"
		return False

	def find_first_level_true_zeroes(self,comp,j,verbose=False):
		bool = True
		zs = self.zeroes(comp,j)
		for z in zs:
			bool = self.find_first_level_true_zero(comp,j,z,verbose=verbose) and bool
		return bool

	def are_zeroes_separated(self,comp,j,verbose=False):
		bool = true
		a = 0
		z = self.zeroes(comp,j)[0]
		while bool and a<len(self.zeroes(comp,j))-1:
			bool = bool and min(self.zero_data[(comp,j,z)][0]) > self.sep_ghost_zero(comp,j,z)
			if verbose and not bool:
				print "Failed at weight",z
			a = a + 1
			z = self.zeroes(comp,j)[a]

		return bool

	def optimize_bounds_naive(self,comp,j):
		p = self.p
		if not self.are_zeroes_separated(comp,j):
			return False
		else:
			zs = self.zeroes(comp,j)
			for z in zs:
				if self.mult(comp,j,z) == 1:
					## ignoring information about total contribution of multiple roots
					d = dimension_cusp_forms(1,z)
					away_from_z = 0
					for zp in zs:
						if zp != z:
							away_from_z += ((ZZ(zp)-ZZ(z)).valuation(p)+1) * self.mult(comp,j,zp)
					self.zero_data[(comp,j,z)][0] = [self.y_lowerbnd(comp,d,z) + (j-d) * (z-2)/2 - away_from_z]
			return True

	def optimize_bounds(self,comp,j):
		p = self.p
		if not self.are_zeroes_separated(comp,j):
			return False
		else:
			zs = self.zeroes(comp,j)
			for z in zs:
				#print "Optimizing at",z
				## ignoring information about total contribution of multiple roots
				d = dimension_cusp_forms(1,z)
				dp = floor(z/3)-1
				away_from_z = 0
				for zp in zs:
					if zp != z:
						away_from_z += ((ZZ(zp)-ZZ(z)).valuation(p)+1) * self.mult(comp,j,zp)
				ptame = []
				## computes valuations coming from level N
				for i in range(d+1):
					if self.good_weight(comp,i,z):
						v = self.y(comp,i,z)
						ptame += [(v,"=")]
					else:
						v = self.y_lowerbnd(comp,i,z)
						ptame += [(v,">=")]
				dnew = dp - 2*d 
				## uses explicit formula for valuations of coefficients arising from p-new forms
				Pnew = Upnew_char_poly(z)
				Pnew_coefs = Pnew.padded_list()
				pnew = [(Pnew_coefs[a].valuation(p),"=") for a in range(Pnew.degree()+1)]
				#print "tame",ptame
				#print "new",pnew
				## forms the product of forms from level N and p-new forms
				assert self.good_weight(comp,d,z), "problem finding classical/non-classical bound"
				bnd = self.y(comp,d,z) + (j-d-1) * (z-2)/2 + (z-1)
				#print "bound",bnd,self.y(comp,d,z),j,d,z
				ans = []
				for a in range(len(ptame)):
					b = j - a
					if b < len(pnew):
						ans += [(ptame[a][0] + pnew[b][0],ptame[a][1])]
						if ans[len(ans)-1][0] > bnd:
							ans[len(ans)-1] = (bnd,">=")
				#print "prod",ans
				ans_vals = [ans[a][0] for a in range(len(ans))]
				m = min(ans_vals)
				if self.mult(comp,j,z) == 1:
					self.zero_data[(comp,j,z)] = [[m - away_from_z], m - away_from_z]
				else:
					self.zero_data[(comp,j,z)][1] = m - away_from_z
			return True

	def find_second_level_true_zero(self,comp,j,z,verbose=False):
		assert self.mult(comp,j,z) > 1,"Not a multiple ghost zero"
		if verbose:
			print "Attempting to find a second zero near",z
		p = self.p
		r = self.sep_nth_level_ghost_zero(comp,j,z,2)
		ks = self.available_sample_weights(comp,j)
		v = [(ZZ(k)-ZZ(z)).valuation(p)+1 for k in ks]
		m = max(v)
		closest = [k for k in ks if (ZZ(k)-ZZ(z)).valuation(p)+1 == m]
		gaps_closest = []
		print "closest",closest
		for k in closest:
			sample = self.prove_yval(comp,j,k) 
			gap = sample[1]
			for zp in self.zeroes(comp,j):
				if not self.close_to_nth_level_ghost_zero(comp,j,k,zp,1):
					gap -= (ZZ(zp)-ZZ(k)).valuation(p) + 1
				else:
					gap = -Infinity
					print "Problem at",k,"with",zp
				if self.mult(comp,j,zp) > 1 and zp != z:
					if self.zero_data[(comp,j,zp)][0][1] > self.sep_nth_level_ghost_zero(comp,j,zp,2):
						gap -= (ZZ(zp)-ZZ(k)).valuation(p) + 1
			gaps_closest += [gap]
		print "gaps",gaps_closest
		e = 1
		print "e",e
		done = false
		while (m-e>=r):
			next_closest = [k for k in ks if (ZZ(k)-ZZ(z)).valuation(p)+1 == m-e]
			print "next_closest",next_closest
			gaps_next_closest = []
			for k in next_closest:
				sample = self.prove_yval(comp,j,k) 
				gap = sample[1]
				for zp in self.zeroes(comp,j):
					if not self.close_to_nth_level_ghost_zero(comp,j,k,zp,1):
						gap -= (ZZ(zp)-ZZ(k)).valuation(p) + 1
					else:
						gap = Infinity
						print "Problem at",k
					if self.mult(comp,j,zp) > 1 and zp != z:
						if self.zero_data[(comp,j,zp)][0][1] > self.sep_nth_level_ghost_zero(comp,j,zp,2):
							gap -= (ZZ(zp)-ZZ(k)).valuation(p) + 1						
				gaps_next_closest += [gap]
			print "next_gaps",gaps_next_closest
			A = max(gaps_closest) 
			B = min(gaps_next_closest)
			if A > B:
				k1 = closest[gaps_closest.index(A)]
				k2 = next_closest[gaps_next_closest.index(B)]
				if verbose:
					print "Gap valuation at",k1,"is",A,"while gap valuation at",k2,"is",B
					print "Success.  There is a true zero closer to",k1,"than to",k2,"\n"
				self.zero_data[(comp,j,z)][0][1] = m-e + .0001
				return True
			else:
				e = e + 1
		print "failed"
		return False
		
	def find_second_level_true_zeroes(self,comp,j,verbose=False):
		bool = True
		zs = self.nth_level_zeroes(comp,j,2)
		for z in zs:
			bool = self.find_second_level_true_zero(comp,j,z,verbose=verbose) and bool
		return bool

	def improve_zero_bounds(self,comp,j,z,verbose=False):
		p = self.p
		r = self.sep_nth_level_ghost_zero(comp,j,z,1)
		ks = self.available_sample_weights(comp,j)
		v = [(ZZ(k)-ZZ(z)).valuation(p)+1 for k in ks]
		m = max(v)
		k = ks[v.index(m)]
		bool = true
		for zp in self.zeroes(comp,j):
			if z != zp:
				bool = bool and (ZZ(k)-ZZ(zp)).valuation(p) + 1 < min(self.zero_data[(comp,j,zp)][0])
			else:
				bnds = self.zero_data[(comp,j,z)][0]
				bool = bool and (ZZ(k)-ZZ(z)).valuation(p) + 1 < min([bnds[a] for a in range(len(bnds)-1)])
		if bool:
			print "Could work!",comp,j,k
		
		sample = self.prove_yval(comp,j,k) 
		gap = sample[1]
		for zp in self.zeroes(comp,j):
			if zp != z:
				gap -= ((ZZ(zp)-ZZ(k)).valuation(p) + 1) * self.mult(comp,j,zp)
			else:
				gap -= ((ZZ(zp)-ZZ(k)).valuation(p) + 1) * (self.mult(comp,j,zp) - 1)
		if gap == m:
			self.zero_data[(comp,j,z)][0][len(self.zero_data[(comp,j,z)][0])-1] = m + .0001
			return True

	def improve_integrality(self,comp,j,z):
		if self.mult(comp,j,z) == 2:
			for a in range(2):
				if floor(self.zero_data[(comp,j,z)][0][a]) != self.zero_data[(comp,j,z)][0][a]:
					if self.zero_data[(comp,j,z)][0][a] * 2 < self.zero_data[(comp,j,z)][1]:
						self.zero_data[(comp,j,z)][0][a] = ceil(self.zero_data[(comp,j,z)][0][a])
		return "Done"

	def gap(self,comp,j,k):
		p = self.p
		sample = self.prove_yval(comp,j,k,verbose=verbose) 
		gap = sample[1]
		for zp in self.zeroes(comp,j):
			if not self.close_to_nth_level_ghost_zero(comp,j,k,zp,1):
				gap -= (ZZ(zp)-ZZ(k)).valuation(p) + 1
			else:
				gap = Infinity
				print "Problem at",k
		return gap


def Upnew_char_poly(k):
	## signs wrong here
	R=PolynomialRing(QQ,'x')
	d = dimension_cusp_forms(1,k)
	dp = floor(k/3)-1
	dnew = dp - 2*d 
	if k%12 == 2 or k%12 == 6:
		return prod([1+x*3**((k-2)/2) for j in range(floor(dnew/2))]) * prod([1-x*3**((k-2)/2) for j in range(ceil(dnew/2))])
	else:
		return prod([1-x*3**((k-2)/2) for j in range(floor(dnew/2))]) * prod([1+x*3**((k-2)/2) for j in range(ceil(dnew/2))])




