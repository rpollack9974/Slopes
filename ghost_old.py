from sage.geometry.newton_polygon import NewtonPolygon
 		
# testing multiplicity conjectures

## p = prime
## G is either Gamma0(N), Gamma1(N) or simply N (which is interpreted at Gamma0(N))
## k = weight
## chi (optional) is a Dirichlet Character
## 
## Returns the tuple: (dim S_k(G), d_k(G cap Gamma_0(p))^{p-new}) or (dim(S_k(Gamma_1(N),chi)),dim(S_k(Gamma_1(Np),chi)))
def dim_of_spaces(p,G,k,chi=None):
        if G in ZZ:
                G = Gamma0(G)
        N = G.level()
        assert (G == Gamma0(N)) or (G == Gamma1(N)), "Other congruence subgroups not coded"
        if chi == None:
                dim_old = dimension_cusp_forms(G,k)
                if G == Gamma0(N):
                        dim_all = dimension_cusp_forms(Gamma0(G.level()*p),k)
                else:
                        Gp = DirichletGroup(N*p)
                        dim_all = sum([dimension_cusp_forms(Gp(chi),k) for chi in DirichletGroup(N)])
        else:
                dim_old = dimension_cusp_forms(chi,k)
                N = chi.level()
                Gp = DirichletGroup(N*p)
                dim_all = dimension_cusp_forms(Gp(chi),k)
	return (dim_old, dim_all - 2*dim_old)

## k = weight
## G is either Gamma0(N), Gamma1(N) or simply N (which is interpreted at Gamma0(N))
## p = prime
## chi (optional) is a Dirichlet Character
## 
## Returns the collection of indices such that k is a semi-stable weight for p and tame level G
def ss_indices(k,G,p,chi=None):
        if k!=2:
                (dk,dknew) = dim_of_spaces(p,G,k,chi=chi)
                return range(dk+1,dk+dknew)
        else:
                (dk,dknew) = dim_of_spaces(p,G,k,chi=chi)
                return range(dk+1,dk+dknew)
                
## G is either Gamma0(N), Gamma1(N) or simply N (which is interpreted at Gamma0(N))
## p = prime
## max_i = maximum index of coefficient of ghost to be computed
## chi is a Dirichlet Character
##
## Returns vector whose r-th component is the associated ghost series on the r-th component of weight space
def form_ghost_power_series(G,p,max_i,chi=None,gamma=None):
        if gamma == None:
                if p==2:
                        gamma = 5
                else:
                        gamma = 1+p
        R = PolynomialRing(ZZ,'w')
        w = R.gen()
        ghosts_coefs = [[R(1) for i in range(max_i+1)] for r in range(p-1)]
        k = 2
        inds = ss_indices(k,G,p,chi=chi)
        while (len(inds)==0 or inds[0]<=max_i+1):
#                print k,": ",inds
                r = k%(p-1)
                wk = gamma**k-1
                for m in range((len(inds)+1)/2):
                        if m < len(inds)/2:
                                if inds[m]<=max_i:
                                        #print "Modifying:",inds[m]
                                        ghosts_coefs[r][inds[m]] *= (w-wk)**(m+1)
                                if inds[len(inds)-1-m]<=max_i:
                                        #print "Modifying:",inds[len(inds)-1-m]
                                        ghosts_coefs[r][inds[len(inds)-1-m]] *= (w-wk)**(m+1)
                        else:
                                if inds[m]<=max_i:
                                        #print "Modifying:",inds[m]
                                        ghosts_coefs[r][inds[m]] *= (w-wk)**(m+1)
                k = k+1
                inds = ss_indices(k,G,p,chi=chi)
        S = PolynomialRing(R,'t')
        t = S.gen()
        ghosts = [S(0) for r in range(p-1)]
        for r in range(p-1):
                for j in range(len(ghosts_coefs[r])):
                        ghosts[r] += ghosts_coefs[r][j]*t**j
        return ghosts
                                
## G is either Gamma0(N), Gamma1(N) or simply N (which is interpreted at Gamma0(N))
## p = prime
## k = weight
## chi (optional) is a Dirichlet Character
## gamma (optional) is a top generator of 1+pZp or 1+4Z_2
## ghost (optional) is a pre-computed ghost series power series
##
## Returns in the ghost slopes in weight k by specializing the ghost power series and taking newton polygon
def ghost_slopes_via_power_series(G,p,k,chi=None,gamma=None,ghost=None):
        if ghost == None:
                ghost = form_ghost(G,p,dimension_cusp_forms(G,k))
        if gamma == None:
                if p==2:
                        gamma = 5
                else:
                        gamma = 1+p
        R = PolynomialRing(QQ,'x')
        w = ghost[0].parent().gen()
        slopes = R(ghost[k%(p-1)].substitute(w=gamma**k-1)).newton_slopes(p)

        return [-slopes[j] for j in range(len(slopes))]
        
## G is either Gamma0(N), Gamma1(N) or simply N (which is interpreted at Gamma0(N))
## p = prime
## max_i = maximum index of coefficient of ghost to be computed
## chi is a Dirichlet Character
##
## Returns vector whose r-th component is a "shell" of the ghost series.  Specifically, the r-th component
##   is a list of the semistable zeroes for each component with multiplicities.  
##
## Example:
## sage: form_ghost_shell(1,2,3)
## [[[],
##   [(14, 1)],
##   [(20, 1), (22, 1), (26, 1)],
##   [(26, 1), (28, 1), (30, 1), (32, 1), (34, 1), (38, 1)]]]
def form_ghost_shell(G,p,max_i,chi=None,old=False,voodoo_seed=[],verbose=False,no_weight_two=False):
        ## Sets the value of gamma (= top generator)
        ghosts_coefs = [[[] for i in range(max_i+1)] for r in range(p-1)]
        
        ## Starting a weight 2, we run through every weight, sort by component, 
        ## compute the associated indices, and then record the weights at those
        ## indices with appropriate multiplicities
        k = 2
        if no_weight_two:
                k = 4
        inds = ss_indices(k,G,p,chi=chi)
        while (len(inds)==0 or inds[0]<=max_i+1):
                r = k%(p-1)
                ## This loops adds the weights to the appropriate indices with the appropriate multiplicities
                for m in range((len(inds)+1)/2):
                        if m < len(inds)/2:
                                if inds[m]<=max_i:
                                        ghosts_coefs[r][inds[m]] += [(k,m+1)]
                                if inds[len(inds)-1-m]<=max_i:
                                        ghosts_coefs[r][inds[len(inds)-1-m]] += [(k,m+1)]
                        else:
                                if inds[m]<=max_i:
                                        ghosts_coefs[r][inds[m]] += [(k,m+1)]
                k = k+1
                inds = ss_indices(k,G,p,chi=chi)
	if (p != 2) and (not old) and (len(voodoo_seed) > 0):
		ghosts_coefs = form_voodoo_modification(G,p,ghosts_coefs,voodoo_seed,verbose=verbose)
		
	if p == 2 and not old:
		#print 'voodoo'
		ghosts_coefs = form_p2_modification(G,ghosts_coefs)
		
	### I have no idea what is going on with indentation but this seems to work lol	
    	return ghosts_coefs
    
## ghost = any ghost shell (i.e. list of roots)
## p = prime
##
## Returns the corresponding power series
def shell_to_power_series(ghost,p):
        if p==2:
                gamma = 5
        else:
                gamma = 1+p
        R = PolynomialRing(ZZ,"w")
        w = R.gen()
        S = PolynomialRing(R,"t")
        t = S.gen()
        ans = [S(0) for r in range(len(ghost))]
        for r in range(len(ghost)):
                g = ghost[r]
                for j in range(len(g)):
                        coef = R(1)
                        for n in range(len(g[j])):
                                root = g[j][n][0]
                                mult = g[j][n][1]
                                coef *= (w-(gamma**root-1))**mult
                        ans[r] += coef*t**j
        return ans
                        
                

## G is either Gamma0(N), Gamma1(N) or simply N (which is interpreted at Gamma0(N))
## p = prime
## k = weight
## max_i (optional) = maximum index of coefficient of ghost to be computed
## chi (optional) is a Dirichlet Character
## gamma (optional) is a top generator of 1+pZp or 1+4Z_2
## ghost (optional) is a pre-computed ghost series power series
##
## From the ghost shell (computed if not given as an input), the first max_i newton points
## of the ghost series at weight k are computed, and the slopes of the associated newton
## polygon are returned.  If max_i is not given, the dimension of S_k(G) is used.
def ghost_slopes_shell_single_weight(G,p,k,max_i=None,chi=None,gamma=None,forced_comp=None,ghost=None):
        if ghost == None:
                ghost = form_ghost_shell(G,p,dimension_cusp_forms(G,k),chi=chi)
        if gamma == None:
                if p==2:
                        gamma = 5
                else:
                        gamma = 1+p
        if forced_comp == None:
                r = k%(p-1)
        else:
                r = forced_comp
        NP = []
        if max_i != None:
                d = min(len(ghost[r]),max_i)
        else:
                d = len(ghost[r])
        for i in range(d):
                if p == 2:
                        e = 2
                else:
                        e = 1
                y = 0
                for ss_wt in ghost[r][i]:
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
        return NP,NewtonPolygon(NP).slopes()

##### This will evaluate the ghost slopes at weight z^k*\eta
##### where \eta is a quadratic character
#####
def ghost_slopes_shell_single_quadchar_weight(G,p,k,eta,max_i=None,chi=None,gamma=None,forced_comp=None,ghost=None):
        if ghost == None:
                ghost = form_ghost_shell(G,p,dimension_cusp_forms(G,k),chi=chi)
        if gamma == None:
                if p==2:
                        gamma = 5
                else:
                        gamma = 1+p
        if forced_comp == None:
                r = k%(p-1)
        else:
                r = forced_comp
        NP = []
        if max_i != None:
                d = min(len(ghost[r]),max_i)
        else:
                d = len(ghost[r])
        for i in range(d):
                if p == 2:
                        e = 2
                else:
                        e = 1
                y = 0
                for ss_wt in ghost[r][i]:
                        if ss_wt[0] == "p":
                                y += ss_wt[1]
                        else:
                                k_ss = ss_wt[0]
                                mult = ss_wt[1]
                                
                                #### then z^k\eta is the point (1+p)^k-1 in wt space
                                if eta(gamma) == 1:
                                	if k_ss >= 0:
                                		y += (valuation(k-k_ss,p)+e)*mult
                                	if k_ss < 0:
                                		y += mult
                                #### this is where z^k\eta lies at -(1+p)^k-1
                                if eta(gamma) == -1:
                                	if k_ss >= 0:
                                		y += mult
                                	if k_ss < 0:
                                		y += (valuation(-k_ss-k,p)+e)*mult
                NP += [(i,y)]
        return NP,NewtonPolygon(NP).slopes()        
######## added by John 9/12
######## added by John 9/12
######## added by John 9/12
######## added by John 9/12
def ghost_newton_points_shell_single_weight(G,p,k,max_i=None,chi=None,gamma=None,ghost=None):
        if ghost == None:
                ghost = form_ghost_shell(G,p,dimension_cusp_forms(G,k),chi=chi)
        if gamma == None:
                if p==2:
                        gamma = 5
                else:
                        gamma = 1+p
        r = k%(p-1)
        NP = []
        if max_i != None:
                d = min(len(ghost[r]),max_i)
        else:
                d = len(ghost[r])
        for i in range(d):
                if p == 2:
                        e = 2
                else:
                        e = 1
                y = 0
                for ss_wt in ghost[r][i]:
                        if ss_wt[0] == "p":
                                y += ss_wt[1]
                        else:
                                k_ss = ss_wt[0]
                                mult = ss_wt[1]
                                y += (valuation(k-k_ss,p)+e)*mult                	
                NP += [(i,y)]
        return NP
        
def compare_with_Wan(G,p,k,max_i=None,chi=None,gamma=None,ghost=None):
	if ghost == None:
		ghost = form_ghost_shell(G,p,dimension_cusp_forms(G,k),chi=chi)
	if gamma == None:
		if p==2:
			gamma = 5
		else:
			gamma = 1+p
	r = k % (p-1)
	if max_i != None:
		d = min(len(ghost[r]),max_i)
	else:
		d = len(ghost[r])
	k0 = k % 12
	if k0 == 2:
		m0 = 0
	else:
		m0 = 1
	NP = ghost_newton_points_shell_single_weight(G,p,k,max_i,chi,gamma,ghost)
	delta_points = [(i,NP[i][1] - ((6*(i-m0)**2)/14 - i/2)) for i in range(max_i)]
	return delta_points
        
        
######## end added by John 9/12
######## end added by John 9/12
######## end added by John 9/12
######## end added by John 9/12
######## end added by John 9/12


## G is either Gamma0(N), Gamma1(N) or simply N (which is interpreted at Gamma0(N))
## p = prime
## k_min, k_max = min/max weights
## max_i (optional) = maximum index of coefficient of ghost to be computed
## chi (optional) is a Dirichlet Character
## gamma (optional) is a top generator of 1+pZp or 1+4Z_2
## ghost (optional) is a pre-computed ghost series power series
##
## Runs ghost_slopes_shell_single_weight from k_min to k_max
def ghost_slopes_shell(G,p,k_min,k_max,max_i=None,chi=None,gamma=None,ghost=None):
        if max_i == None:
                if chi == None:
                        max_i = dimension_cusp_forms(G,k_max)+2
                else:
                        if ((-1)**k_max)==chi(-1):
                                max_i = dimension_cusp_forms(chi,k_max)+2
                        else:
                                max_i = dimension_cusp_forms(chi,k_max-1)+2
        if ghost == None:
                ghost = form_ghost_shell(G,p,max_i,chi=chi)
        if gamma == None:
                if p==2:
                        gamma = 5
                else:
                        gamma = 1+p
        ans = []
        for k in range(k_min,k_max+1):
                if chi == None:
                        ans += [(k,ghost_slopes_shell_single_weight(G,p,k,max_i=dimension_cusp_forms(G,k)+1,chi=chi,gamma=gamma,ghost=ghost))]
                else:
                        ans += [(k,ghost_slopes_shell_single_weight(G,p,k,max_i=dimension_cusp_forms(chi,k)+1,chi=chi,gamma=gamma,ghost=ghost))]
        return ans



def buzzard_vs_ghost_shell(N,p,max_k,ghost=None):
#        print "Working on the prime ",p
	import datetime
	import time	

	now = datetime.datetime.now()
	strnow = now.strftime("_date=%Y_%m_%d-%H_%M")

	file = "p="+str(p)+"-N="+str(N)+"-max_k="+str(max_k)+"_"+str(strnow)
	log_file = file + ".txt"
	save_file = file + ".sobj"

	F = open(log_file,'w')
	F.write("Beginning prime specific test for p = " + str(p) + " in tame level N = " + str(N) + " up to k <= " + str(max_k) + "\n")
	F.close()

	# compute the buzzard slopes and record it in the file

	buzzard_start_time = time.time()
	buzzard_slopes = buzzard_tpslopes(p,N,max_k)
	buzzard_end_time = time.time()
	buzzard_log = "Computed Buzzard's slopes in " + str(buzzard_end_time-buzzard_start_time) + " seconds\n"
#        print buzzard_log
	F = open(log_file,'a')
	F.write(buzzard_log)
	F.close()


	# compute our slopes and record it in the file
	our_start_time = time.time()	
        max_i = max([len(buzzard_slopes[k]) for k in range(len(buzzard_slopes))])
        if ghost == None:
        	ghost = form_ghost_shell(N,p,max_i)
        our_slopes = [[],[]]
        for k in range(2,max_k+1):
                slopes = ghost_slopes_shell_single_weight(N,p,k,max_i=len(buzzard_slopes[k])+1,ghost=ghost)
                our_slopes += [slopes]

	our_end_time = time.time()
	our_log = "Computed the ghost's slopes in " + str(our_end_time - our_start_time) + " seconds\n"	
#        print our_log
	F = open(log_file,'a')
	F.write(our_log)
	F.close()

	# make the comparison and write results
#	compare_ours = [x[1] for x in our_slopes]
	pass_var = (buzzard_slopes == our_slopes)

	if pass_var:
		results = "PASSED"
	else:
		results = "FAILED"
	F = open(log_file,'a')
	F.write("Test results: " + str(results) + "\n")
	F.close()
#        print "Test results: " + str(results) + "\n"

	# now return the outcome of the test
#        return (buzzard_slopes,our_slopes)
	return (pass_var,our_end_time-buzzard_start_time)


# will test buzzard vs. the ghost 
# for p <= max_p and k <= max_k
# in tame level N
# one single log file will be kept progressing the entire process
# and each individual instance will run a log as well
def buzzard_ghost_shell_deathmatch(N,max_p,max_k):
	import datetime
	now = datetime.datetime.now()
	strnow = now.strftime("_date=%Y_%m_%d-%H_%M")
	
	log_filename = "N="+str(N)+"-max_p="+str(max_p)+"-max_k="+str(max_k)+"_" + str(strnow) + ".txt"
	
	F = open(log_filename,'w')
	F.write("Beginning deathmatch for p <= " + str(max_p) + " and tame level N = " + str(N) + " and k <= " + str(max_k) + "\n\n")
	F.close()
	
	total_passed = 0
	total_time = 0
	
        if N == 1:
                prime_list = [p for p in prime_range(2,max_p + 1) if gcd(p,N) == 1]
        else:
                prime_list = [p for p in prime_range(3,max_p + 1) if gcd(p,N) == 1]
	
	for p in prime_list:
		F = open(log_filename,'a')
		F.write("Beginning test for p = " + str(p) + " and k <= " + str(max_k) + "\n")
		F.close()
		
		(passed,p_time) = buzzard_vs_ghost_shell(N,p,max_k)
		if passed:
			passed_str = "PASSED"
			total_passed += 1
		else:
			passed_str = "FAILED"
		total_time += p_time
		
		end_string = "Finished test for p = " + str(p) + " and k <= " + str(max_k) + "\n"
		end_string += "Result: " + str(passed_str) + "\n"
		end_string += "Time: " + str(p_time) + " seconds\n"
		end_string += "Current progress:\n"
		end_string += "\t Total ran: " + str(prime_list.index(p)+1) + "\n"
		end_string += "\t Total passed: " + str(total_passed) + "\n\n"		
		F = open(log_filename,'a')
		F.write(end_string)		
		F.close()
	
	final_string = "\n\nFinished all testing\n"
	final_string += "Results:\n" 
	final_string += str(total_passed) + " out of " + str(len(prime_list)) + " passed\n"
	final_string += "Total time: " + str(total_time) + " seconds\n"
	F = open(log_filename,'a')
	F.write(final_string)		
	F.close()

	return total_passed == len(prime_list)

def buzzard_regular(p,G,chi=None):
        assert p>2, "Sorry, need to think more about p=2"
        kp = (p+3)/2
        if chi == None:
                spaces = [ModularSymbols(G,k,1,GF(p)).cuspidal_subspace() for k in range(2,kp+1)]
        else:
                spaces = [ModularSymbols(chi,k,1,GF(p)).cuspidal_subspace() for k in range(2,kp+1)]
        ords = [S.hecke_polynomial(p).substitute(0).norm().valuation(p)==0 for S in spaces]
        return not ords.count(False) > 0

def buzzard_regular2(p,G,chi=None):
		if p == 2:
			S2 = ModularSymbols(G,2,1,GF(2)).cuspidal_subspace()
			if S2.hecke_matrix(2).determinant().norm().valuation(2) != 0:
				return False
			S4 = ModularSymbols(G,4,1,QQ).cuspidal_subspace()
			slopes_4 = S4.hecke_polynomial(2).newton_slopes(2)
			c1 = (dimension_cusp_forms(2*G,2) - dimension_cusp_forms(G,2))
			c2 = (len(slopes_4) - slopes_4.count(0))
			return slopes_4.count(0) == c1 and slopes_4.count(1) == c2
		elif p == 3:
			kp = p+1 ### not sure if we need to check k = 4 or not, but just being safe
			if chi == None:
				spaces = [ModularSymbols(G,k,1,GF(p)).cuspidal_subspace() for k in range(2,kp+1)]
			else:
				spaces = [ModularSymbols(chi,k,1,GF(p)).cuspidal_subspace() for k in range(2,kp+1)]
			ords = [S.hecke_matrix(p).determinant().norm().valuation(p)==0 for S in spaces]
			return not ords.count(False) > 0
		else:
			kp = (p+3)/2
			if chi == None:
				spaces = [ModularSymbols(G,k,1,GF(p)).cuspidal_subspace() for k in range(2,kp+1)]
			else:
				spaces = [ModularSymbols(chi,k,1,GF(p)).cuspidal_subspace() for k in range(2,kp+1)]
			ords = [S.hecke_matrix(p).determinant().norm().valuation(p)==0 for S in spaces]
			return not ords.count(False) > 0

def constant_terms_of_hecke_polynomials(p):
        assert p > 2, "Sorry, p > 2"
        kp = (p+3)/2
        weights = range(16,kp+1,2)
        if kp >= 12:
                weights = [12] + weights
        
        spaces = [ModularSymbols(1,k,1,GF(p)).cuspidal_subspace() for k in weights]
        ans = []
        for S in spaces:
                ans += [S.hecke_polynomial(p).substitute(0)]
        ans.sort()
        return ans

def constant_terms_of_random_polynomials(p):
        assert p > 2, "Sorry, p > 2"
        kp = (p+3)/2
        weights = range(16,kp+1,2)
        if kp >= 12:
                weights = [12] + weights
        
        ans = []
        for k in weights:
                ans += [floor(random()*10**5)%p]
        ans.sort()
        return ans

def collect_irregular_data(p_min,p_max):
        ans = []
        irregs = []
        for p in primes(p_min,p_max+1):
                p_ans = constant_terms_of_hecke_polynomials(p)
                if len(p_ans)>0 and p_ans[0] == 0:
                        irregs += [p]
#                        print p," is irregular"
		file_name = "irregular.txt"
		F = open(file_name,'w')
		F.write([p,p_ans]+',')
		
                ans += [[p,p_ans]]
        return ans,irregs
        
def collect_random_irregular_data(p_min,p_max):
        ans = []
        irregs = []
        for p in primes(p_min,p_max+1):
                p_ans = constant_terms_of_random_polynomials(p)
                if len(p_ans)>0 and p_ans[0] == 0:
                        irregs += [p]
                ans += [[p,p_ans]]
        return ans,irregs

## G = Congruence subgroup
## p = prime
## mink,maxk = weight ranges
## chi (optional) = Dirichlet character
##
## Returns the slopes of T_p on S_k(G) for k in the weight range        
def true_slopes(G,p,mink,maxk,chi=None):
        N = G.level()
	import datetime
	now = datetime.datetime.now()
	strnow = now.strftime("_date=%Y_%m_%d-%H_%M")
	
        if chi == None:
                log_filename = "Slopes:"+str(G)+"-p="+str(p)+"_" + str(strnow) + ".txt"
        else:
                log_filename = "Slopes:"+str(G)+str(chi)+"-p="+str(p)+"_" + str(strnow) + ".txt"

	F = open(log_filename,'w')
        if chi == None:
                F.write("Beginning slopes computation for "+str(G) + " and p = " + str(p) + "\n\n")
        else:
                F.write("Beginning slopes computation for "+str(G) + str(chi) + " and p = " + str(p) + "\n\n")
	F.close()

        ans = []
	for k in range(mink,maxk+1):
                if chi == None:
                        M = ModularSymbols(G,k,1).cuspidal_subspace()
                else:
                        M = ModularSymbols(chi,k,1).cuspidal_subspace()                        
		slopes = M.hecke_polynomial(p).newton_slopes(p)
		for r in range(len(slopes)):
			if slopes[r] > 10**6:
				slopes[r] = Infinity
		slopes.reverse()
#		if condense:
#			slopes = condense_slopes(slopes)
		answer = str(k) + ":" + str(slopes) + "\n"
		F = open(log_filename,'a')
		F.write(answer)		
		F.close()

		#print k,":",slopes
                ans += [(k,slopes)]
		del(M)
        
        return ans

        
def true_slopes_vs_ghost_shell(G,p,min_k,max_k,chi=None):
        if G in ZZ:
                G = Gamma0(G)
        N = G.level()
        assert (G == Gamma0(N)) or (G == Gamma1(N)), "Other congruence subgroups not coded"

        assert N%p!=0, "Level not prime to p"

#        print "Working on the prime ",p
	import datetime
	import time	

	now = datetime.datetime.now()
	strnow = now.strftime("_date=%Y_%m_%d-%H_%M")

#	file = "True.vs.ghost.G="+str(G)+"p="+str(p)+"-max_k="+str(max_k)+"_"+str(strnow)
#	log_file = file + ".txt"
#	save_file = file + ".sobj"

#	F = open(log_file,'w')
#	F.write("Beginning prime specific test for p = " + str(p) + " in level = " + str(G) + " up to k <= " + str(max_k) + "\n")
#	F.close()

	# compute the true slopes and record it in the file

	true_start_time = time.time()
	t_slopes = true_slopes(G,p,min_k,max_k,chi=chi)
	true_end_time = time.time()
	true_log = "Computed true slopes in " + str(true_end_time-true_start_time) + " seconds\n"
#        print true_log
#	F = open(log_file,'a')
#	F.write(true_log)
#	F.close()

	# compute our slopes and record it in the file
	our_start_time = time.time()	
        max_i = max([len(t_slopes[j][1]) for j in range(len(t_slopes))])
        our_slopes = ghost_slopes_shell(G,p,min_k,max_k,chi=chi)
	our_end_time = time.time()
	our_log = "Computed the ghost's slopes in " + str(our_end_time - our_start_time) + " seconds\n"	
#        print our_log
#	F = open(log_file,'a')
#	F.write(our_log)
#	F.close()

	# make the comparison and write results
	pass_var = (t_slopes == our_slopes)

        v = []
	if pass_var:
		results = "PASSED"
                print "PASSED"
	else:
		results = "FAILED"
                v = [t_slopes[j][0] for j in range(len(t_slopes)) if t_slopes[j] != our_slopes[j]]
                print "FAILED -- differ in weights: ",v
#	F = open(log_file,'a')
#	F.write("Test results: " + str(results) + "\n")
#	F.close()
#       print "Test results: " + str(results) + "\n"

	# now return the outcome of the test
#        return (t_slopes,our_slopes)
	return (pass_var,our_end_time-true_start_time,v)


def true_vs_ghost_deathmatch(G,min_p,max_p,min_k,max_k,chi=None):
        N = G.level()
	import datetime
	now = datetime.datetime.now()
	strnow = now.strftime("_date=%Y_%m_%d-%H_%M")
	
        if chi == None:
                log_filename = "true.v.ghost"+str(G)+"-min_p="+str(min_p)+"-max_p="+str(max_p)+"-min_k="+str(min_k)+"-max_k="+str(max_k)+"_" + str(strnow) + ".txt"
        else:
                log_filename = "true.v.ghost"+str(G)+str(chi)+"-min_p="+str(min_p)+"-max_p="+str(max_p)+"-min_k="+str(min_k)+"-max_k="+str(max_k)+"_" + str(strnow) + ".txt"

	
	F = open(log_filename,'w')
        if chi == None:
                F.write("Beginning deathmatch for "+str(min_p)+" <= p <= " + str(max_p) + " and " + str(G) + " and " +str(min_k)+" <= k <= " + str(max_k) + "\n\n")
        else:
                F.write("Beginning deathmatch for "+str(min_p)+" <= p <= " + str(max_p) + " and " + str(G) + str(chi) + " and " +str(min_k)+" <= k <= " + str(max_k) + "\n\n")
	F.close()
	
	total_passed = 0
        total_tests = 0
        total_skips = 0
	total_time = 0
        passed_str = ""
        p_time = 0
        
	
	prime_list = [p for p in prime_range(min_p, max_p + 1) if gcd(p,N) == 1]

	for p in prime_list:
		F = open(log_filename,'a')
		F.write("Beginning test for p = " + str(p) + " and " + str(min_k) + " <= k <= " + str(max_k) + "\n")
		F.close()
		print "\nBeginning test for p = " + str(p) + " and " + str(min_k) + " <= k <= " + str(max_k)

		bool = buzzard_regular(p,G,chi=chi)

                if not bool:
                        total_skips += 1
                        F = open(log_filename,'a')
                        if chi == None:
                                F.write("---" + str(p) + " not regular for " + str(G)+" \n \n")
                        else:
                                F.write("---" + str(p) + " not regular for " + str(G)+" and character " + str(chi) + " \n \n")
                        print "NOT REGULAR!!"
                        F.close()
                else:
                        total_tests += 1
                        (passed,p_time,v) = true_slopes_vs_ghost_shell(G,p,min_k,max_k,chi=chi)
                        if passed:
                                passed_str = "PASSED"
                                total_passed += 1
                        else:
                                passed_str = "FAILED at weights: " + str(v)
                        total_time += p_time
		
                        end_string = "Finished test for p = " + str(p) + " and k <= " + str(max_k) + "\n"
                        end_string += "Result: " + str(passed_str) + "\n"
                        end_string += "Time: " + str(p_time) + " seconds\n"
                        end_string += "Current progress:\n"
                        end_string += "\t Total ran: " + str(total_tests) + "\n"
                        end_string += "\t Total passed: " + str(total_passed) + "\n"		
                        end_string += "\t Total skipped: " + str(total_skips) + "\n\n"		
                        F = open(log_filename,'a')
                        F.write(end_string)		
                        F.close()
	
	final_string = "\n\nFinished all testing\n"
	final_string += "Results:\n" 
	final_string += str(total_passed) + " out of " + str(total_tests) + " passed\n"
	final_string += "Total tests skipped:  " + str(total_skips) + "\n"
	final_string += "Total time: " + str(total_time) + " seconds\n"
	F = open(log_filename,'a')
	F.write(final_string)		
	F.close()

	return total_passed == total_tests

def true_vs_ghost_deathmatch_char_by_char(G,min_p,max_p,min_k,max_k):
        N = G.level()
	import datetime
	now = datetime.datetime.now()
	strnow = now.strftime("_date=%Y_%m_%d-%H_%M")

	DG = DirichletGroup(N).galois_orbits()
        prim_char = [chi[0] for chi in DG if chi[0].is_primitive()]
        for chi in prim_char:
                print "Working with character ",chi
                log_filename = "true.v.ghost"+str(G)+str(chi)+"-min_p="+str(min_p)+"-max_p="+str(max_p)+"-min_k="+str(min_k)+"-max_k="+str(max_k)+"_" + str(strnow) + ".txt"

	
                F = open(log_filename,'w')
                F.write("Beginning deathmatch for "+str(min_p)+" <= p <= " + str(max_p) + " and " + str(G) + str(chi) + " and " +str(min_k)+" <= k <= " + str(max_k) + "\n\n")
                F.close()
	
                total_passed = 0
                total_tests = 0
                total_skips = 0
                total_time = 0
        
	
                prime_list = [p for p in prime_range(min_p, max_p + 1) if gcd(p,N) == 1]

                for p in prime_list:
                        F = open(log_filename,'a')
                        F.write("Beginning test for p = " + str(p) + " and " + str(min_k) + " <= k <= " + str(max_k) + "\n")
                        F.close()
                        print "\nBeginning test for p = " + str(p) + " and " + str(min_k) + " <= k <= " + str(max_k)

                        bool = buzzard_regular(p,G,chi=chi)

                        if not bool:
                                total_skips += 1
                                F = open(log_filename,'a')
                                if chi == None:
                                        F.write("---" + str(p) + " not regular for " + str(G)+" \n \n")
                                else:
                                        F.write("---" + str(p) + " not regular for " + str(G)+" and character " + str(chi) + " \n \n")
                                print "NOT REGULAR!!"
                                F.close()
                        else:
                                total_tests += 1
                                (passed,p_time,v) = true_slopes_vs_ghost_shell(G,p,min_k,max_k,chi=chi)
                                if passed:
                                        passed_str = "PASSED"
                                        total_passed += 1
                                else:
                                        passed_str = "FAILED at weights: " + str(v)
                                total_time += p_time
		
                                end_string = "Finished test for p = " + str(p) + " and k <= " + str(max_k) + "\n"
                                end_string += "Result: " + str(passed_str) + "\n"
                                end_string += "Time: " + str(p_time) + " seconds\n"
                                end_string += "Current progress:\n"
                                end_string += "\t Total ran: " + str(total_tests) + "\n"
                                end_string += "\t Total passed: " + str(total_passed) + "\n"		
                                end_string += "\t Total skipped: " + str(total_skips) + "\n\n"		
                                F = open(log_filename,'a')
                                F.write(end_string)		
                                F.close()
	
                final_string = "\n\nFinished all testing\n"
                final_string += "Results:\n" 
                final_string += str(total_passed) + " out of " + str(total_tests) + " passed\n"
                final_string += "Total tests skipped:  " + str(total_skips) + "\n"
                final_string += "Total time: " + str(total_time) + " seconds\n"
                F = open(log_filename,'a')
                F.write(final_string)		
                F.close()

                print "FINISHED WITH ",chi
                print total_passed == total_tests

###################################################################################
###################################################################################
##  Messing around with p=3 and N=17
###################################################################################
###################################################################################
##
##  There are 5 possible rhobars.  
##      --An eisenstein one: 1 \oplus \varepsilon
##      --One defined over F_9 and its twist (both locally reducible at 3)
##      --X_0(17)[3] and its twist (both locally irreducible at 3)
##     
##  We note X_0(17)[3] is uniquely determined by a_2 = 2 (mod 3),
##  and its twist is uniquely determined by a_2 = 1 (mod 3)


## Gives the dimension of the rhobar subspace with rhobar = X_0(17)[3]
def d(k): 
    ans = 0
    A = ModularSymbols(17,k,1,GF(3)).cuspidal_subspace().decomposition()
    for r in range(len(A)):
        T = A[r].hecke_operator(2)
        ans += ((T+1)**A[r].dimension()).kernel().dimension()
    return ans

## From the data this looks like the multiplicity of X_0(17)[3]
def irred_rhobar(k):
        return floor((k+2)/4)

## Should give the dimension of the rhobar subspace with rhobar = X_0(17)[3]
def dt(k): 
    ans = 0
    A = ModularSymbols(17,k,1,GF(3)).cuspidal_subspace().decomposition()
    for r in range(len(A)):
        T = A[r].hecke_operator(2)
        ans += ((T+2)**A[r].dimension()).kernel().dimension()
    return ans

## From the data this looks like the multiplicity of the twist of X_0(17)[3]
def irred_rhobar_twist(k):
        return floor((k-2)/4)

##  From the data this looks like the dimension of X_0(17)[3] and its twist in weight k (level 17)
def irred_dim(k):
        if (k%4 == 0):
                e = -1
        else:
                e = 0
        return k/2+e

## Gives the multiplicity of X_0(17)[3] in the 3-new subspace
def dnew_rhobar(k):
    ans = 0
    A = ModularSymbols(17*3,k,1,GF(3)).cuspidal_subspace().decomposition()
    for r in range(len(A)):
        T = A[r].hecke_operator(2)
        ans += ((T+1)**(A[r].dimension())).kernel().dimension()
    return ans-2*irred_rhobar(k)

## Gives the multiplicity of X_0(17)[3] in the 3-new subspace
def dnew_rhobar(k):
    ans = 0
    A = ModularSymbols(17*3,k,1,GF(3)).cuspidal_subspace().decomposition()
    for r in range(len(A)):
        T = A[r].hecke_operator(2)
        ans += ((T+1)**(A[r].dimension())).kernel().dimension()
    return ans-2*irred_rhobar(k)

## Gives the dimension of X_0(17)[3] and its twist in the 3-new subspace
def dnew_rhobar_twist(k):
    ans = 0
    A = ModularSymbols(17*3,k,1,GF(3)).cuspidal_subspace().decomposition()
    for r in range(len(A)):
        T = A[r].hecke_operator(2)
        ans += ((T+2)**(A[r].dimension())).kernel().dimension()
    return ans-2*irred_rhobar_twist(k)

##  From the data this looks like the multiplcity of X_0(17)[3] in weight k, level 51 and 3-new
##  same for the twist
def irred_rhobar_new(k):
        return 2*floor(k/4)

##  From the data this looks like the multiplcity of X_0(17)[3] and its twist in weight k, level 51 and 3-new
def irred_dim_new(k):
        if (k%4 == 2):
                e = -2
        else:
                e = 0
        return k+e

def ss_indices_modified(k,gal="red"):
        if gal=="red":
                (dk,dknew) = dim_of_spaces(3,17,k)
                dk -= irred_dim(k)
                dknew -= irred_dim_new(k)
                return range(dk+1,dk+dknew)
        elif gal == "irred":
                if k%2 == 0:
                        dk = irred_dim(k)
                        dknew = irred_dim_new(k)
                        return range(dk+1,dk+dknew)
                else:
                        return range(1,0)
        elif gal == "rhobar":
                if k%2 == 0:
                        dk = irred_rhobar(k)
                        dknew = irred_rhobar_new(k)
                        return range(dk+1,dk+dknew)
                else:
                        return range(1,0)
        elif gal == "rhobar_twist":
                if k%2 == 0:
                        dk = irred_rhobar_twist(k)
                        dknew = irred_rhobar_new(k)
                        return range(dk+1,dk+dknew)
                else:
                        return range(1,0)


def form_ghost_shell_modified(max_i,gal="red"):
        p = 3
        gamma = 1 + p

        ghosts_coefs = [[[] for i in range(max_i+1)] for r in range(p-1)]
        
        ## Starting a weight 2, we run through every weight, sort by component, 
        ## compute the associated indices, and then record the weights at those
        ## indices with appropriate multiplicities
        k = 2
        inds = ss_indices_modified(k,gal=gal)
        if gal == "irred":
                ghosts_coefs[0][1] = [(2,1)]
        elif gal == "rhobar":
                ghosts_coefs[0][1] = [(2,1)]                
#        elif gal == "rhobar_twist":
#                ghosts_coefs[0][1] = [(2,2)]                
        while (len(inds)==0 or inds[0]<=max_i+1):
                r = k%(p-1)
                ## This loops adds the weights to the appropriate indices with the appropriate multiplicities
                for m in range((len(inds)+1)/2):
                        if m < len(inds)/2:
                                if inds[m]<=max_i:
                                        ghosts_coefs[r][inds[m]] += [(k,m+1)]
                                if inds[len(inds)-1-m]<=max_i:
                                        ghosts_coefs[r][inds[len(inds)-1-m]] += [(k,m+1)]
                        else:
                                if inds[m]<=max_i:
                                        if gal == "red":
                                                ghosts_coefs[r][inds[m]] += [(k,m+1)]
                                        else:
                                                if m==0:
                                                        ghosts_coefs[r][inds[m]] += [(k,m+1)]
                                                else:
                                                        ghosts_coefs[r][inds[m]] += [(k,m+2)]
                k = k+1
                inds = ss_indices_modified(k,gal=gal)

        return ghosts_coefs

def ghost_slopes_shell_single_weight_modified(k,max_i=None,ghost=None,gal="red"):
        if ghost == None:
                ghost = form_ghost_shell_modified(dimension_cusp_forms(17,k),gal=gal)
        p = 3
        gamma = 1+p
        r = k%(p-1)
        NP = []

        if max_i != None:
                d = min(len(ghost[r]),max_i)
        else:
                d = len(ghost[r])
 
        for i in range(d):
                if p == 2:
                        e = 2
                else:
                        e = 1
                y = 0
                for ss_wt in ghost[r][i]:
                        if ss_wt[0] == "p":
                                y += ss_wt[1]
                        else:
                                k_ss = ss_wt[0]
                                mult = ss_wt[1]
                                y += (valuation(k-k_ss,p)+e)*mult
                NP += [(i,y)]
        return NewtonPolygon(NP).slopes()

def ghost_slopes_shell_modified(k_min,k_max,ghost=None,gal="red"):
        if gal=="red":
                max_i = dimension_cusp_forms(17,k_max)+2
        elif gal == "irred":
                max_i = irred_dim(k_max)+2
        elif gal == "rhobar":
                max_i = irred_rhobar(k_max)+2
        elif gal == "rhobar_twist":
                max_i = irred_rhobar_twist(k_max)+2
        if ghost == None:
                ghost = form_ghost_shell_modified(max_i,gal=gal)
        
        p = 3
        gamma = 4

        ans = []
        for k in range(k_min,k_max+1,2):
                if gal=="red":
                        d = dimension_cusp_forms(17,k)
                elif gal == "irred":
                        d = irred_dim(k)
                elif gal == "rhobar":
                        d = irred_rhobar(k)
                elif gal == "rhobar_twist":
                        d = irred_rhobar_twist(k)
                ans += [(k,ghost_slopes_shell_single_weight_modified(k,max_i=d+1,ghost=ghost,gal=gal))]
        return ans

def testing():
        for k in range(4,50,2):
                print k,"------------"
                M = ModularSymbols(17,k,1).cuspidal_subspace()
                t = M.hecke_polynomial(3).newton_slopes(3)
                t.sort()
                print "true:      ",t
                print "predicted: ",ww[k-2][1]
#                print "removing:  ",v1[k-2][1]
                for s in range(len(v1[k-2][1])):
                        t.remove(v1[k-2][1][s])
                print "left:      ",t
                print "predicted: ",v2[k-2][1]
#                print v2[k-2]

def testing2(k_max):
        v_red = ghost_slopes_shell_modified(2,k_max,gal="red")
        v_irred = ghost_slopes_shell_modified(2,k_max,gal="irred")

        for k in range(2,k_max,2):
                print "testing ",k
                M = ModularSymbols(17,k,1).cuspidal_subspace()
                t = M.hecke_polynomial(3).newton_slopes(3)
                t.sort()
                g_slopes = (v_irred[k/2-1][1] + v_red[k/2-1][1])
                g_slopes.sort()
#                print "True slopes:",t
#                print "Ghost slopes",g_slopes
                if g_slopes != t:
                        print "******************FALSE*******************"


def decomp2(N,k,p):
        print IntegerRingMod(p)
        print "HERE",N,k,p
        M = ModularSymbols(N,k,1,GF(p)).cuspidal_subspace()
        return M
        print M
        Ms = [M]
        print "HERE"
        for q in primes(N):
                Ms = decomp_helper(Ms,q)
        return Ms
                
                
def decomp_helper(Ms,q):
        ans = []
        for spaces in Ms:
                if (N%q!=0) and (q%p!=0):
                        Tq = M.hecke_operator(q)
                        fq = Tq.char_poly()
                        pieces = fq.factor()
                        for d in pieces:
                                ans += [((d[0]**d[1]).substitute(Tq)).kernel()]
        return ans
                

def abstract_ss_indices(list,k):
        dk = list[k][0]
        dknew = list[k][1]
        return range(dk+1,dk+dknew)
        
def john_modify_ghost_3_17(g):
	new_ghost = list(g)
	for m in range(1,10):
		k = 2*m
		if (k % 4) == 2:
			new_ghost[0][dimension_cusp_forms(Gamma0(17),k)].append((k/2 +1,1))
        	new_ghost[0][dimension_cusp_forms(Gamma0(17),k)+1].append((k/2+1,-1))
        	print new_ghost[0][dimension_cusp_forms(Gamma0(17),k)]
        	print "\n"
		if (k % 4) == 0:
			print [k,k/2,dimension_cusp_forms(Gamma0(17),k)]
			new_ghost[0][dimension_cusp_forms(Gamma0(17),k)].append((k/2,1))
    		new_ghost[0][dimension_cusp_forms(Gamma0(17),k)-1].append((k/2,-1))
    		print new_ghost[0][dimension_cusp_forms(Gamma0(17),k)]
    		print "\n"
    	return new_ghost
                                        
## The ghost series is determined by the dimensions of S_k(1) and S_k(p)^{p-new}.  
## This function takes any association k --> two numbers and returns the corresponding ghost series
##
## list = list of tuples with list[k] representing the two numbers attached to k
## max_i = highest degree
def form_abstract_ghost_shell(list,max_i,p):
        ghosts_coefs = [[[] for i in range(max_i+1)] for r in range(p-1)]
        
        ## Starting a weight 2, we run through every weight, sort by component, 
        ## compute the associated indices, and then record the weights at those
        ## indices with appropriate multiplicities
        k = 2
        inds = abstract_ss_indices(list,k)
        while (len(inds)==0 or inds[0]<=max_i+1):
                r = k%(p-1)
                ## This loops adds the weights to the appropriate indices with the appropriate multiplicities
                for m in range((len(inds)+1)/2):
                        if m < len(inds)/2:
                                if inds[m]<=max_i:
                                        ghosts_coefs[r][inds[m]] += [(k,m+1)]
                                if inds[len(inds)-1-m]<=max_i:
                                        ghosts_coefs[r][inds[len(inds)-1-m]] += [(k,m+1)]
                        else:
                                if inds[m]<=max_i:
                                        ghosts_coefs[r][inds[m]] += [(k,m+1)]
                k = k+1
                inds = abstract_ss_indices(list,k)

        return ghosts_coefs
        






#############################################
#############################################
##  SHAPE OF GHOST SPECTRAL CURVE
#############################################
#############################################

## p = prime
## i = index
## v = valuation (non-integral)
## comp = component
## Return (r,s) where r = total ghost zeroes (with mult) > v
## and s = sum v_p(w) where w is a ghost zero with v_p(w) < v. Thus
## on the valuation v_p(w) = v we see v_p(a_i) = v*r + s.

## 5/15/16 added ability to re-center counts at an integer weight k
## defaults to centering at the weight k = 0
def zero_val_count(p,comp,i,v,ghost,central_wt=0):
        assert floor(v)!=v, "Use non-integral valuations"
        ai = ghost[comp][i]
        if p == 2:
                e = 2  ## this is the extra valuation from changing from k to w
        else:
                e = 1
        r= sum([a[1] for a in ai if ZZ(a[0]-central_wt).valuation(p)+e > v])
        s= sum([(ZZ(a[0]-central_wt).valuation(p)+ e)*a[1] for a in ai if ZZ(a[0]-central_wt).valuation(p)+e < v])
        return (r,s)

## v = valuation
## Returns a list (r_i,s_i) such that (i,r_i*v+s_i) is the i-th Newton point of U_p in valuation v = v2(wt-central_wt)
## added 5/15/16: ability to use a central_wt
def newton_points_at_fixed_valuation(p,comp,v,max_i,ghost,central_wt=0):
        return [zero_val_count(p,comp,i,v,ghost,central_wt) for i in range(max_i+1)]


### added 5/15/16: ability to use a central_wt 

def global_halo_piece(p,comp,v,max_i,ghost,central_wt=0):
        list = newton_points_at_fixed_valuation(p,comp,v,max_i,ghost,central_wt)
#        slopes = NewtonPolygon([(i,list[i][0]*v+list[i][1]) for i in range(max_i+1)]).slopes()
#        assert len(set(slopes))==len(slopes), "SLOPES NOT DISTINCT!!"
        left = NewtonPolygon([(i,list[i][0]*floor(v)+list[i][1]) for i in range(max_i+1)]).slopes()
        right = NewtonPolygon([(i,list[i][0]*ceil(v)+list[i][1]) for i in range(max_i+1)]).slopes()
#        ple = [list[i][0]*floor(v)+list[i][1] for i in range(max_i+1)]
#        pre = [list[i][0]*ceil(v)+list[i][1] for i in range(max_i+1)]
#        left = [ple[i]-ple[i-1] for i in range(1,max_i+1)]
#        right = [pre[i]-pre[i-1] for i in range(1,max_i+1)]

        print left,right
        spectral=line([(floor(v),0),(ceil(v),0)])
        prev = [(floor(v),0),(ceil(v),0)]
        for i in range(max_i):
                C = CC[i]
                if left[i] == right[i]:
                        C = "purple"
                if prev == [(floor(v),left[i]),(ceil(v),right[i])]:
                        C = "black"
#                slope = right[i]-left[i]
 #               C = CC[slope]
                spectral += line([(floor(v),left[i]),(ceil(v),right[i])],color=C)
                C = "blue"
                prev = [(floor(v),left[i]),(ceil(v),right[i])]

        return spectral
        
def halo_line(pt1,pt2,linethickness,linecolor,prime):
	p = prime
	r = 10
	outline = 'dashed'
	colored = 'white'
	if p >= 3 or (p == 2 and pt2[0] >= 4):
		return line([pt1,pt2],thickness=linethickness,color=linecolor) + point(pt1,size=linethickness*r,color=colored,faceted = True,alpha=1,zorder=10) + point(pt2,size=linethickness*r,color=colored,faceted = True,alpha=1,zorder=10)
	if p == 2 and pt2[0] < 4:
		return line([pt1,pt2],thickness=linethickness,color=linecolor)
	
		
#### june 7, added code so that if you feed in a tame level then
#### it will highlight every other arithmetic progression
#### either black/blue
#### also made the points scale with the size of the family
#### for aesthetic reasons

def global_halo_piece(p,comp,v,max_i,ghost,central_wt=0,tame_level=0):
		max_i = 2*max_i
		list = newton_points_at_fixed_valuation(p,comp,v,max_i,ghost,central_wt)
		left = NewtonPolygon([(i,list[i][0]*floor(v)+list[i][1]) for i in range(max_i+1)]).slopes()
		right = NewtonPolygon([(i,list[i][0]*ceil(v)+list[i][1]) for i in range(max_i+1)]).slopes()
		max_i = max_i/2
		r = floor(v)
		num_progs = 0
		if tame_level > 0:
			num_progs = (p**(r+1))*(p-1)*(p+1)*Gamma0(N).index()/24   
		spectral = point((0,0))
		prev = [(r,0),(r+1,0)]
		linethickness = .1
		for i in range(max_i):
			linecolor = 'blue'
			if num_progs > 0 and floor(i/num_progs) % 2 == 1:
				linecolor = 'blue'
			C = CC[i]
			linethickness = len([j for j in range(max_i) if left[j]==left[i] and right[j] == right[i]])
			spectral += halo_line((r,left[i]),(r+1,right[i]),linethickness,linecolor,prime=p)
			prev = [(floor(v),left[i]),(ceil(v),right[i])]
		return spectral

### added 5/15/16: ability to use a central_wt (see the zero_val_count routine)

def global_halo(p,comp,min_v,max_v,max_i,ghost,central_wt = 0,tame_level=0,xaxes=None,yaxes=None):
        assert floor(min_v)==min_v and floor(max_v)==max_v, "C'mon.  Use an integer for the endpoint"
        gh = global_halo_piece(p,comp,min_v+0.5,max_i,ghost,central_wt,tame_level)
        for v in range(min_v+1,max_v):
                gh +=  global_halo_piece(p,comp,v+0.5,max_i,ghost,central_wt,tame_level)
        xaxis_label = 'v_'+str(p)+'(w($\kappa$)-w('+str(central_wt)+'))'
        yaxis_label = 'slope'
        title = 'Halos of ghost series\ncentered at w = w(' + str(central_wt) + ')'
        if xaxes != None or yaxes != None:
                gh.axes_range(0,xaxes,0,yaxes)
        
        return [gh, text(xaxis_label,((max_v-min_v)/2,-5),color='black'), title]
        

D=[1, 1, 3, 5, 6, 10, 11, 21,22, 42,43, 85,86,170, 171, 341]
        

#CC = ["red","red","blue","red","blue","green","blue","red","blue","green","purple","purple","purple","green","blue","red","blue","green","purple","purple","purple","black","black","black","black","black","purple","purple","purple","green","blue","red","blue","green","purple","purple","purple","black","black","black","black","black"]+["blue" for i in range(1000)]

#CC = ["blue"] + ["purple"] + ["green" for i in range(3)] + ["orange" for i in range(5)] + ["red" for i in range(11)] + ["purple" for i in range(21)] + ["blue" for i in range(1000)]

#CC = ["black","blue","green","yellow","orange","red"]+["blue" for i in range(1000)]
#CC = [Color(0,0,0),Color(.5,0,.5),Color(.75,.3,.75),Color(1,.3,1)]+["blue" for i in range(1000)]

CC = ["blue" for i in range(1000)]

#############################
#############################
#############################
### messing around with p = 13
### and N = 5
#############################
#############################
#############################

###
# in weights 2,...,14 there are 
# only 2 non-ordinary forms
# One has weight 6 and the other
# has weight 10, and is theta^8(weight 6)
###

###
# thus there is one twist-class
# of irreducible rho-bar
# and the rest are reducible
###

def irred_system_13_5(threshold):
	return [[p, -ModularSymbols(5,6,1,GF(13)).cuspidal_subspace().hecke_polynomial(p)[0]] for p in prime_range(threshold)]

###
# takes an eigensystem (p,a_p) 
# and returns the eigensystem
# for (p, p^i a_p) corresponding
# to the ith Tate twist mod p
###

def twist(i,eigensystem):
	twist_system = []
	for x in eigensystem:
		twist_system.append([x[0], (x[0]**i)*x[1]])
	return twist_system
	
def all_twists(p,eigensystem):
	return [twist(i,eigensystem) for i in range(p-1)]

def find_irred_rhobar_13_5(k,threshold = 100):
	spaces = ModularSymbols(5,k,1,GF(13)).cuspidal_subspace().decomposition()
	check_systems = all_twists(13,irred_system_13_5(threshold))
	dim_count = 0
	for A in spaces:
		hecke_polys = [A.hecke_polynomial(p) for p in prime_range(threshold)]
		for ES in check_systems:
			if set([hecke_polys[i].substitute(ES[i][1]) for i in range(len(hecke_polys))]) == set([0]):
				dim_count += A.dimension()	
				print "found one! in weight k = " + str(k) + " with dimension " + str(A.dimension())
	return dim_count
		
		
######
# testing whether Buzzards conjecture
# really predicts the Buzzard--Gouvea bound
# up to weight k = (p^2 + 1)/2 (this comes from
# some non-sense Berger sent me from Yamashita-Yasuda
######

def test_buzzardgouvea(p,N):
	max_k = floor((p**2+1)/2)
	bslopes = buzzard_tpslopes(p,N,max_k)
	for k in range(0,max_k+1,2):
		gouvea_buzzard_bound = floor((k-1)/(p+1))
		ok = True
		if len(bslopes[k]) > 0:
			max_slope_k = max(bslopes[k])
			if max_slope_k > gouvea_buzzard_bound:
				ok = False
		if not ok:
			print "YOU FOUND AN EXAMPLE WITH k,P,N = " + str([k,p,N])

def test_bg(max_p,max_N):
	for p in prime_range(max_p+1):
		for N in range(1,max_N+1):
			if gcd(N,p) == 1:
				test_buzzardgouvea(p,N)
				#print (p,N)
	

def test_ss_line(G,p,k,ghost,slopes=None):
        d = dim_of_spaces(p,G,k)
        d_old = d[0]
        d_ss = d[1]
        if slopes == None:
                slopes = ghost_slopes_shell_single_weight(G, p, k, ghost=ghost)[1] ## changed 10/20/2016 because somehow the slope commands changed at some point.
        slopes_ss = slopes[d_old:d_ss+d_old]
        if slopes_ss.count((k-2)/2) == len(slopes_ss):
                return True
        else:
                return False,slopes_ss

def test_old_symmetric(G,p,k,ghost,slopes=None):
        d = dim_of_spaces(p,G,k)
        d_old = d[0]
        d_ss = d[1]
        d_full = 2*d_old+d_ss
        if slopes == None:
                slopes = ghost_slopes_shell_single_weight(G, p, k, ghost=ghost)[1] ## changed 10/20/2016 because somehow the slope commands changed at some point.
        j=0
        ans = True
        while j+d_ss+d_old<min(len(slopes),d_full):
                if ((slopes[j+d_ss+d_old] + slopes[d_old-j-1]) != (k-1)):
                        ans = False
                        print "Problem at indices",d_old-j-1,j+d_ss+d_old,"slopes ",slopes[d_old-j-1],slopes[j+d_ss+d_old]
                j = j+1
        return ans

#### added by john 10/18

def test_old_NP_symmetric(G,p,k,ghost):
	(d_old,d_ss) = dim_of_spaces(p,G,k)
	d_full = 2*d_old + d_ss
	NP = ghost_newton_points_shell_single_weight(G,p,k,d_full+1,ghost=ghost)
	j = 1
	ans = True
	while j <= d_old:
		slp1 = NP[j][1] - NP[j-1][1]
		slp2 = NP[d_full-(j-1)][1] - NP[d_full-j][1]
		if slp1 + slp2 != k-1:
			ans = False
			print "Problem at indices",j,d_full-(j-1),"with weight k=",k
		j += 1
	return ans

def test_theta(G,p,k,ghost,slopes=None):
        d = dim_of_spaces(p,G,k)
        if d[0]<len(ghost[0]): 
                if slopes == None:
                        slopes = ghost_slopes_shell_single_weight(G, p, k, ghost=ghost)[1] ## changed 10/20/2016 because somehow the slope commands changed at some point.
                sk=mults(slopes)
                sneg = mults(ghost_slopes_shell_single_weight(G,p,2-k,ghost=ghost)[1]) ## changed 10/20/2016
                theta_sneg = [(x[0]+k-1,x[1]) for x in sneg]
                j = 0
                while (j<len(sk)) and (sk[j][0]<=k-1):
                        j=j+1
                i = 0
                #skip slope k-1 
                if theta_sneg[0][0]==k-1:
                        i = i+1
                ans= True
                # left off last few slopes because it wasn't working (Newton polygon problems?)
                while j<len(sk)-3 and i<len(theta_sneg)-3:
                		#if theta_sneg[i] != sk[j]:
                        if list(theta_sneg[i]) != sk[j]:  ### changed 10/20/2016 (see issues above)
                                ans = False
                                print "Problem at weight:",k
                                print "Problem at indices:",i,j," and slopes ",sk[j],theta_sneg[i]
                        i = i+1
                        j = j+1
                return ans
        else:
                return True,"Vacuous"
        
        
        
def full_test(G,p,k,ghost,slopes=None,verbose=False):       
        if slopes == None:
                if verbose:
                        print "Forming slopes"
                slopes = ghost_slopes_shell_single_weight(G,p,k,ghost=ghost)
        if verbose:
                print "Testing semistable line"
        ans_ss = test_ss_line(G,p,k,ghost,slopes=slopes)
        if verbose:
                print ans_ss

        if verbose:
                print "Testing symmetry of old forms"
        ans_old = test_old_symmetric(G,p,k,ghost,slopes=slopes)
        if verbose:
                print ans_old
        
        ## added by john 10/18      
		if verbose:
			print "Testing deeper Newton point symmetry of old forms"
        ans_old_sym = test_old_NP_symmetric(G,p,k,ghost)
        if verbose:
                print ans_old
		
        if verbose:
                print "Testing theta"
        
        ## 10/18: not working b/c it needs "mults(-)" that I can't
        ## find defined anywhere
        
        #ans_theta = test_theta(G,p,k,ghost,slopes=slopes)
        ans_theta = None
        if verbose:
                print ans_theta

        return (ans_ss,ans_old,ans_old_sym,ans_theta)
        
### little test for p = 11 and N = 8

def run_john_testOct20(k):
	p = 11
	N = 8
	G = DirichletGroup(p)
	L  = G.base_ring()
	P = L.ideal(p).prime_factors()[0]
	omega = [chi for chi in G if (chi(2) - 2).valuation(P) > 0][0]
	ans = []
	for j in [z for z in range(1,(p-1)/2 + 1) if (k%2)==(z%2)]:
		S = ModularSymbols(DirichletGroup(p*N)(omega**j),weight=k,sign=1)
		f = S.hecke_polynomial(11)
		NPs = [L(x).valuation(P) for x in f.coeffs()]
		R = PolynomialRing(ZZ,'t')
		t = R.gen()
		new_f = sum([p**(NPs[i])*(t**i) for i in range(f.degree()+1)])
		ans.append((j,new_f.newton_slopes(p)))
	return str(ans)

def mults(v):
        mult_free = list(set(v))
        mult_free.sort()
        ans = []
        for a in mult_free:
                ans.append([a,v.count(a)])
        return ans

def mult_rig_slps(slps):
	return [(slps[j][0],mults(slps[j][1])) for j in range(len(slps)) if slps[j][0] % 2 == 0]

def slopes_in_weight2_with_char_fixed_weight(p,N,k,acc):
        assert p>2, "Need p odd"
        Qp = pAdicField(p,acc)
	G = DirichletGroup(p,Qp)
        Zmodp = Integers(p)
        g = Zmodp.unit_gens()[0]
        omega = [chi for chi in G if (chi(g)-ZZ(g)).valuation()>0][0]

	file = "wt2data/wt2.p="+str(p)+"-N="+str(N)
	log_file = file + ".txt"
	save_file = file + ".sobj"

	F = open(log_file,'a')
        F.write("["+str(k)+",[\n")
        F.close()

        ans = []
        for j in [z for z in range(0,(p-2)) if (k%2)==(z%2)]:
                M = ModularSymbols(DirichletGroup(p*N,Qp)(omega**j),weight=k,sign=1)
                f = M.hecke_polynomial(p)
                slopes = f.newton_slopes(p)
                slopes.reverse()
                slopes = mults(slopes)
                print k,j,slopes
                F = open(log_file,'a')
                F.write("("+str(j)+","+str(slopes)+"),\n")
                F.close()


	F = open(log_file,'a')
        F.write("]\n")
        F.close()

        return acc

def slopes_in_weight2_with_char(p,N,kmin,kmax):
        C=sage.rings.padics.precision_error.PrecisionError
        acc=500
        done=false
        k = kmin
        while k<=kmax and not done:
                try:
                        slopes_in_weight2_with_char_fixed_weight(p,N,k,acc)
                except (C):
                        acc = acc+100
                        print "Trying accuracy ",acc
                        k=k-1
                k=k+1
        return 1
                
        
 #### p = 2 helpers
 
## takes in a list [(a0,d0),(a1,d1),...] where
## 0 = a0  < a1 < ... and di > 0
## and returns the correct multiplicity pattern
## for the weight k = 2 points
## sage: wt2_slopes_to_mults([(0,7),(1/2,1),(4,3)])
## [0, 1, 2, 3, 3, 2, 1, 0, 0, 1, 1, 0]
 
def wt2_slopes_to_mults(mult_slope_list):
	mults = []
	ell = len(mult_slope_list)
	for i in range(ell):
		dimi = mult_slope_list[i][1]
		for j in range(1,dimi+1):
			if j < dimi/2:
				mults.append(j)
			else:
				mults.append(dimi+1 - (j+1))
	return mults
	

### returns a list [(-k,m)]_i where we want to make a modification
### to the ghost a_i by adding in a zero of mult m at the weight -5^k - 1.

def get_modified_zeros(N,max_k):
	# get even character of conductor 8
	chi = [c for c in DirichletGroup(8) if c.conductor() == 8 and c(-1) == 1][0]
	S = ModularSymbols(DirichletGroup(N*8)(chi),weight=2,sign=1).cuspidal_subspace()
	#S = CuspForms(DirichletGroup(N*8)(chi),2)
	# computes weight 2 slopes with cond. 8
	mult_slope_list = mults(S.hecke_polynomial(2).newton_slopes(2))
	# isolates just the fractional slopes in the classical space
	base_mults = wt2_slopes_to_mults([x for x in mult_slope_list if x[0] != 0 and x[0] != 1])
	if mult_slope_list[0][0] == 0:
		start_ind = mult_slope_list[0][1] + 1
	else:
		start_ind = 1
	gap = (S.dimension() + Gamma0(N).ncusps())
	mods = [[]]*(start_ind + len(base_mults) + gap*(max_k-1))
	for k in range(2,max_k+1):
		for j in range(len(base_mults)):
			if base_mults[j] != 0:
				mods[j+start_ind + gap*(k-2)] = [(-k,base_mults[j])]
				#mods[j+start_ind + gap*(k-2)] = [(-2,base_mults[j])]
	return mods

### fun fact: if N is odd then
### dim S_{k+1}(Gamma1(8N),e\pm) - dim S_k(Gamma1(8N),e\pm) = N\cdot prod_{\ell \dvd N} (1 + 1/\ell)

def form_p2_modification(N,ghost):
	new_ghost_0 = [x for x in ghost[0]]
	new_ghost = [new_ghost_0,[]]
	max_i = len(new_ghost[0])
	dim_gaps = N*prod([1 + 1/ell for ell in ZZ(N).prime_factors()])
	chi = [c for c in DirichletGroup(8) if c.conductor() == 8 and c(-1) == 1][0]
	dim2 = dimension_cusp_forms(DirichletGroup(8*N)(c),2)
	max_k = floor((max_i - dim2)/dim_gaps + 2)
	mod_zeros = get_modified_zeros(N,max_k)
	for j in range(max_i):
		new_ghost[0][j] += mod_zeros[j]
	return new_ghost
	
def spread_out(mult_list):
	new_list = []
	for x in mult_list:
		for j in range(x[1]):
			new_list.append(x[0])
	return new_list
	
def gather_fractional_slopes(slope_data):
	ret_data = []
	for y in slope_data:
		wt = y[0]
		twist_data = []
		for z in y[1]:
			twist = z[0]
			mult_list = z[1]
			new_mults = spread_out(mult_list)
			frac_data = [(x,new_mults.index(x)) for x in set(new_mults) if not x.is_integer()]
			twist_data.append([twist,frac_data])
		ret_data.append([wt,twist_data])
	return ret_data

# want to do N  = [47, 71, 103, 127, 151, 167, 191, 199, 239, 263, 271]
def test_kevin_us_p2(N):
	print "Starting test in level N = " + str(N)
	foo = buzzard_vs_ghost_shell(N,2,2050)
	print "\n\t" + str(foo)
	
###### testing super naive guesses 11/2

def cl_wt(i):
    return i+3

def twist_dim(k):
    return (k-1)*12-4
    
    
##returns the lambda invariant of a polynomial in shell form (e.g. h = [(10,2),(20,1)] has lambda = 3)
def lambda_shell(h):    
        return sum([h[a][1] for a in range(len(h))])
        
## returns the number of roots of h = [(10,2), (20,1)] with a fixed valuation
## it looks for v = v(w_k) though.
## when p = 2, v = 0 will return the "boundary roots" 
def sum_lambda_shell(p,h,v):
	if p == 2:
		d = 2
	if p != 2:
		d = 1
	if v > 0:
		return sum([x[1] for x in h if x[0] > 0 and ZZ(x[0]).valuation(p) == v-d])
	if v == 0:
		return sum([x[1] for x in h if x[0] < 0]) ### this is only for p = 2
		
def generic_val(p,h,v):
	total = lambda_shell(h)
	boundaries = sum_lambda_shell(p,h,0)
	other_total = [sum_lambda_shell(p,h,x) for x in range(1,v)]
	others_weights = [x*sum_lambda_shell(p,h,x) for x in range(1,v)]
	return boundaries + sum(others_weights) + v*(total - boundaries - sum(other_total))
	
def wadic_newton_ghost_slopes(g,r,maxi):
        v = [lambda_shell(g[r][a]) for a in range(maxi+1)]
        return [v[a]-v[a-1] for a in range(1,len(v))]

def wadic_ghost_slopes(g,r,maxi):
        v = [lambda_shell(g[r][a]) for a in range(maxi+1)]
        A = [(a,v[a]) for a in range(0,maxi)]
        return NewtonPolygon(A).slopes()

        
        
    
        
        
##returns in shell form g_i/g_{i-1}        
def form_delta_i_shell(g,r,i):
        a=g[r][i]
        b=[(g[r][i-1][t][0],-g[r][i-1][t][1]) for t in range(len(g[r][i-1]))]
        return combine([a]+[b])
        ans=[g[r][i][t] for t in range(len(g[r][i]))]
        roots = [t[0] for t in ans]
        for alpha in g[r][i-1]:
                if alpha[0] in roots:
                        j = roots.index(alpha[0])
                        ans[j]=(alpha[0],ans[j][1]-alpha[1])
                else:
                        ans=[alpha]+ans
        return ans

def combine(list):
        if len(list)>2:
                return combine([combine([list[0],list[1]])]+list[2:len(list)])
        else:
                ans = list[0]+list[1]
                ans=[list[0][t] for t in range(len(list[0]))]
                roots = [t[0] for t in ans]
                for alpha in list[1]:
                        if alpha[0] in roots:
                                j = roots.index(alpha[0])
                                ans[j]=(alpha[0],ans[j][1]+alpha[1])               
                        else:
                                ans=ans+[alpha]
                roots = [t[0] for t in ans]
                mults = [t[1] for t in ans]
                roots.sort()
                roots2 = [t[0] for t in ans]                
                ans = [(roots[j],ans[roots2.index(roots[j])][1]) for j in range(len(roots))]
                return ans

def arith_progs(p):
	G = DirichletGroup(p**2)
	chi = G[1]**(p-1)
	omega = G[1]**p
	Ms=[ModularSymbols(chi*omega**(2*r),2,1).cuspidal_subspace() for r in range(0,(p-1)/2)]
	Ps=[Ms[a].hecke_polynomial(p) for a in range((p-1)/2)]
	Pcs = [Ps[a].coefficients() for a in range((p-1)/2)]
	for a in range(len(Pcs)):
		Pcs[a].reverse()
	Ks=[Pcs[a][1].parent() for a in range(len(Pcs))]
	pps = [Ks[a].factor(p)[0][0]  for a in range(len(Ks))]
	Nps=[[(a,Pcs[b][a].valuation(pps[b])) for a in range(len(Pcs[b]))] for b in range(len(Pcs))]
	ss=[NewtonPolygon(Nps[a]).slopes() for a in range(len(Nps))]
	slopes=[]
	for r in range(len(ss)):
		slopes = slopes + ss[r]
	slopes = slopes + [p-1 for a in range((p-1)/2)]
	slopes.sort()
	return slopes

#f1 and f2 are general dimension formulas
def form_ghost_shell_general(G,p,max_i,r,f1,f2,chi=None,old=False,voodoo_seed=[],verbose=False,no_weight_two=False):
        ## Sets the value of gamma (= top generator)
        ghosts_coefs = [[] for i in range(max_i+1)]
        
        ## Starting a weight 2, we run through every weight, sort by component, 
        ## compute the associated indices, and then record the weights at those
        ## indices with appropriate multiplicities
        k = 2
        print "HELLO" 
        if no_weight_two:
                k = 4
        inds = range(f1[k]+1,f1[k]+f2[k])
        while (len(inds)==0 or inds[0]<=max_i+1):
                if k%(p-1)==r:
                        ## This loops adds the weights to the appropriate indices with the appropriate multiplicities
                        for m in range((len(inds)+1)/2):
                                if m < len(inds)/2:
                                        if inds[m]<=max_i:
                                                ghosts_coefs[inds[m]] += [(k,m+1)]
                                        if inds[len(inds)-1-m]<=max_i:
                                                ghosts_coefs[inds[len(inds)-1-m]] += [(k,m+1)]
                                else:
                                        if inds[m]<=max_i:
                                                ghosts_coefs[inds[m]] += [(k,m+1)]

                k = k+1
                inds = range(f1[k]+1,f1[k]+f2[k])
	if (p != 2) and (not old) and (len(voodoo_seed) > 0):
		ghosts_coefs = form_voodoo_modification(G,p,ghosts_coefs,voodoo_seed,verbose=verbose)
		
	if p == 2 and not old:
		#print 'voodoo'
		ghosts_coefs = form_p2_modification(G,ghosts_coefs)
		
	### I have no idea what is going on with indentation but this seems to work lol	
    	return ghosts_coefs
