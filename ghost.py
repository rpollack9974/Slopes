from sage.geometry.newton_polygon import NewtonPolygon

class ghost(SageObject):
	"""This class represents the ghost series for a prime p and a fixed tame level Gamma_0(N).
		By ghost series, we mean the (p-1)/2 power series for each even component of weight space

		To initialize a ghost series simply type:

			ghost(p,N) 

		or to specify the number of coefficients computed:

			ghost(p,N,num_coefs=??)
		
		The key functions of this class are:

		1) slopes
		2) wadic_slopes

		If g is the ghost series, the first is called by:

			g.slopes(k) or g.slopes(k,num=??) 

		and returns the ghost slopes in weight k.

		The second is called by:

			g.wadic_slopes(comp) or g.wadic_slopes(comp,num=??)

		and returns the w-adic slopes on the component comp.
	"""

	def compute_ghost_series(self,num_coefs=10):
		"""This function computes and stores the relevant ghost seris in 'shell' form --- that is,
			for each component comp, it sets self.series[comp] to a list whose i-th term is a 
			list of the ghost zeroes of the i-th coefficient with multiplicities.

			Inputs:
				num_coefs (optional) -- an integer representing how many coefficients are computed.
		"""
		p=self.p
		N=self.N

		ghost_coefs = [[[] for i in range(num_coefs+1)] for a in range(p-1)]

		## Starting at weight 2, we run through weights in the component,
		## compute the associated indices, and then record the weights at those
		## indices with appropriate multiplicities

		k=2;
		if k==0:
			k=k+p-1
		f1 = dimension_cusp_forms(N,k)
		f2 = dimension_cusp_forms(N*p,k)-2*f1

		inds = range(f1+1,f1+f2)
		while (len(inds)==0 or inds[0]<=num_coefs+1):
			## This loops adds the weights to the appropriate indices with the appropriate multiplicities
			for m in range((len(inds)+1)/2):
				if m < floor(len(inds)/2):
					if inds[m]<=num_coefs:
						ghost_coefs[k%(p-1)][inds[m]] += [(k,m+1)]
					if (inds[len(inds)-1-m]<=num_coefs):
						ghost_coefs[k%(p-1)][inds[len(inds)-1-m]] += [(k,m+1)]
				else:
					if inds[m]<=num_coefs:
						ghost_coefs[k%(p-1)][inds[m]] += [(k,m+1)]
			k = k+1
			f1 = dimension_cusp_forms(N,k)
			f2 = dimension_cusp_forms(N*p,k)-2*f1
			inds = range(f1+1,f1+f2)
		self.series=ghost_coefs

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

	def __repr__(self):
		return "Ghost series for p="+str(self.p)+" and N="+str(self.N)

	def __getitem__(self,a):
		"""Returns the power series (in shell form) on a-th component"""
		return self.series[a]

	def ghost_zeroes(self,comp,i):
		"""Returns a list of the ghost series of the i-th coefficient on component comp"""
		return self.series[comp][i]

	def num_coefs(self,comp):
		"""Returns the number of coeffcients computed on component comp"""
		return len(self[comp])

	def slopes(self,k,num=None):
		"""Returns the slopes in weight k --- num-many or all computed if term==None

		HACKING IN THIS ONE ON HOW MANY TERMS TO USE"""
		NP = []
		p=self.p
		comp=k%(p-1)
		if num==None:
			d = self.num_coefs(comp)
		else:
			### HACKING HERE (UNCLEAR HOW MANY TERMS NEED TO BE USED TO GET CORRECT SLOPES)
			if self.num_coefs(comp)<num+10:
				self.compute_ghost_series(num_coefs=num+10)
			d = min(self.num_coefs(comp),num+10) 
		if p == 2:
			e = 2
		else:
			e = 1
		for i in range(d):
			y = 0
			for ss_wt in self[comp][i]:
				if ss_wt[0] == "p":
					y += ss_wt[1]
				else:
					k_ss = ss_wt[0]
					mult = ss_wt[1]
						
				#### added by john 10/17, see form_ghost_shell for instructions
				if k_ss >= 0:
					y += (valuation(k-k_ss,p)+e)*mult
				if k_ss < 0:
					y += mult
			NP += [(i,y)]

		if num==None:
			return NewtonPolygon(NP).slopes()
		else:
			return NewtonPolygon(NP).slopes()[0:num]

	def multiplicity(self,comp,i):
		"""Returns the total number of zeroes with multiplicity in the i-th coefficient
			on component comp"""
		return sum([self[comp][i][a][1] for a in range(len(self[comp][i]))])

	def wadic_slopes(self,comp,num=None):
		"""Returns the w-adic slopes of the mod p reduction of the ghost series on 
			component comp"""
		if num!=None:
			##HACKING HERE ABOUT NUMBER OF TERMS
			if self.num_coefs(comp)<num+10:
				g.compute_ghost_series(num_coefs=num+10)

		NP = [(a,self.multiplicity(comp,a)) for a in range(self.num_coefs(comp))]
		if num!=None:
			return NewtonPolygon(NP).slopes()[0:num]
		else:
			return NewtonPolygon(NP).slopes()

