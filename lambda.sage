def dims_on_rim(p):
	G = DirichletGroup(p^2)
	psi = (G.0)^(p-1)
	omega = (G.0)^(p)
	return [ModularSymbols(psi*omega^(2*r-2),2,1).cuspidal_subspace().dimension() for r in range(1,(p-1)/2+1)]
	
def slopes_on_rim(N,p,comp):
	assert p%2==1, "Need an odd prime.  Sorry."
	G = DirichletGroup(p^2)
	GN = DirichletGroup(N*p^2)
	psi = G.0
	chi = psi^(p-1)
	K = psi.base_ring()
	pp = K.prime_above(p)
	e = pp.ramification_index()
	omega = psi^p
	r=0
	v = [(omega(a)^r-a).valuation(pp) for a in range(1,p)]
	while min(v)==0:
		r=r+1
		v = [(omega(a)^r-a).valuation(pp) for a in range(1,p)]
	omega = omega^r
	assert omega.order() == p-1
	slopes = [a*(p-1) for a in range(1,(p+1)/2)] #eisenstein stuff
	for j in range(1,(p+1)/2):
#		print "Computing the ",j,"-th modular symbol space of ",GN((chi*omega^(comp-2))^(-1)*omega^(2*(j-1)))
		M = ModularSymbols(GN((chi*omega^(comp-2))^(-1)*omega^(2*(j-1))),2,1,base_ring=K).cuspidal_subspace()
#		print "---computing Hecke"
		f = M.hecke_polynomial(p)
		f = f.coefficients()
		f.reverse()
		v = [(r,f[r].valuation(pp)/e) for r in range(len(f))]
		v = NewtonPolygon(v).slopes()
#		print v
		v = [(j-v[a]) for a in range(len(v))]
#		print "shifted",v
		v = [v[a]*(p-1) for a in range(len(v))]
#		print "scaled",v 
		slopes = slopes + v
		slopes.sort()
#		print "Slopes so far:",slopes
	
	return slopes

def mult(g,j,comp):
    return sum([g[comp][j][i][1] for i in range(len(g[comp][j]))])


def test_on_rim(p,comp):
    g=form_ghost_shell(1,p,p*(p^2-1)/24+30)
    v1 = NewtonPolygon([(j,mult(g,j,comp)) for j in range(0,p*(p^2-1)/24+1)]).slopes()
    v2 = slopes_on_rim(1,p,comp)
    return v1 == v2
	


########## Trying to find mod p representaitons (level N)
def find_mod_q_eigenforms(p,N):
    ans = []
    for k in range(2,p^2+2):
    	M=ModularSymbols(N,k,1,GF(p))
	M = M.decomposition()
	for r in range(len(M)):
	    Tqs=[]
	    for q in primes(20):
	    	if (N*p)%q != 0:
	    	   Tqs += [(q,M[r].hecke_matrix(q)[0,0])]
	    	ans += [Tqs]
    return ans
	 

def cut_down(fs):
    v=[]
    for r in range(len(fs)):
    	if v.count(fs[r])==0:
	   v+=[fs[r]]
    return v
	    

def counting_fixed_eigensystem_fixed_weight(p,N,k,aps):
    M=ModularSymbols(N,k,1,GF(p)).cuspidal_subspace()
    for r in range(len(aps)):
    	q=aps[r][0]
	Tq=M.hecke_operator(q)
	M=((Tq-aps[r][1])^10).kernel()
    return M.dimension()

def counting_fixed_eigensystem(p,N,kmin,kmax,aps,comp=-1):
#    v=[0,0]
    v=[]
    for k in range(kmin,kmax,1):
    	if comp==-1 or (k%(p-1))==comp:
    	   print k
    	   v+=[counting_fixed_eigensystem_fixed_weight(p,N,k,aps)]
	   print "---",v[len(v)-1]
    return v

#def counting_fixed_new_eigensystem_fixed_weight(p,N,k,aps):
#    M=ModularSymbols(p*N,k,1,GF(p)).cuspidal_subspace().new_subspace()
#    for r in range(len(aps)):
#    	q=aps[r][0]
#	Tq=M.hecke_operator(q)
#	M=((Tq-aps[r][1])^10).kernel()
#    return M.dimension()

#def counting_fixed_new_eigensystem(p,N,kmin,kmax,aps):
#    v=[0,0]
#    for k in range(kmin,kmax,1):
#    	v+=[counting_fixed_new_eigensystem_fixed_weight(p,N,k,aps)]
#    return v
	   
def counting_fixed_levelp_eigensystem_fixed_weight(p,N,k,aps):
    M=ModularSymbols(p*N,k,1,GF(p)).cuspidal_subspace()
    for r in range(len(aps)):
    	q=aps[r][0]
	Tq=M.hecke_operator(q)
	M=((Tq-aps[r][1])^10).kernel()
    return M.dimension()

def counting_fixed_levelp_eigensystem(p,N,kmin,kmax,aps,comp=-1):
#    v=[0,0]
    v=[]
    for k in range(kmin,kmax,1):
    	if comp==-1 or (k%(p-1))==comp:
	   print k
    	   v+=[(k,counting_fixed_levelp_eigensystem_fixed_weight(p,N,k,aps))]
	   print "---",v[len(v)-1]
    return v

#def counting_fixed_levelp_eigensystem(p,N,kmin,kmax,aps,comp=-1):
#    v=[0,0]
#    for k in range(kmin,kmax,1):
#    	if comp!=-1:
#	   if k%(p-1)==comp:
#	       	print k
#    		v+=[counting_fixed_levelp_eigensystem_fixed_weight(p,N,k,aps)]
#		for j in range(p-2):
#		    v+=[0]
#	else:
#	   print k
#    	   v+=[counting_fixed_levelp_eigensystem_fixed_weight(p,N,k,aps)]
#    return v
	   
#p=11, N=1 v is regular dims, w is full dim (not new)	
#sage: v=[]
#sage: for k in range(0,10000):
#    if k%10==2:
#        v+=[floor((k-2)/120)+1]
#    else:
#        v+=[0]

#sage: w=[]
#sage: for k in range(0,10000):
#    if k%10==2:
#        w+=[floor((k-2)/10)-1] 
#    else:
#        w+=[0]
#....:      


def level1_irred(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[floor((k-12)/120)+1]
	else:
	   v+=[0]
    return v

def new_irred(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[floor((k-2)/10)+1-2*(floor((k-12)/120)+1)]
	else:
	   v+=[0]
    return v

def full_irred(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[floor((k-2)/10)+1]
	else:
	   v+=[0]
    return v

def level1_irred_twist(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[floor((k-2+50)/120)]
	else:
	   v+=[0]
    return v

def new_irred_twist(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[ceil((k-2)/10)-2*floor((k-2+50)/120)]
	else:
	   v+=[0]
    return v

def full_irred_twist(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[ceil((k-2)/10)]
	else:
	   v+=[0]
    return v

def level1_1_10(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[floor((k-12)/120)]
	else:
	   v+=[0]
    return v

def new_1_10(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[floor((k-22)/10)+1-2*(floor((k-12)/120))]
	else:
	   v+=[0]
    return v

def full_1_10(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[floor((k-22)/10)+1]
	else:
	   v+=[0]
    return v

def level1_2_9(N):
    return augment(11,seqn(11,8,4,3,floor(N/(p-1))),2)



def full_2_9(N):
    v=[]
    for k in range(0,N):
    	if k%10==2:
	   v+=[(k-2)/5]
	else:
	   v+=[0]
    return v    

def new_2_9(N):
    a=level1_2_9(N)
    b=full_2_9(N)
    return [b[r]-2*a[r] for r in range(min(len(a),len(b)))]

def seqn(p,a,b,s,N):
    v=[0 for i in range(s)]
    for j in range(a):
    	v+=[1]
    for j in range(b):
    	v+=[2]
    i=a+b+s
    while i<N:
    	  v+=[v[len(v)-(p+1)]+2]
	  i+=1
    return v

def augment(p,seq,r):
    v=[0 for j in range(len(seq)*(p-1)+r)]
    for j in range(len(seq)):
    	v[r+j*(p-1)]=seq[j]
    return v
    
def ab(p,a,b,s,r,N):
    return augment(p,seqn(p,a,b,s,N),r)

def abp(p,a,b,s,r,N):
    v=seqn(p,a,b,s,N)
    w=[2*(i+1)-2*v[i] for i in range(N)]
    return augment(p,w,r)

def level1_13_delta(N):
    augment(13,seqn(13,12,2,1,floor(N/(p-1)),0))

def level32_5(N):
    v=[1,1,3,3,3,3]
    for j in range(N/len(v)):
    	for i in range(6):
	    v+=[v[len(v)-6]+4]
    return v

def level5_32_5(N):
    v=[2+4*i for i in range(N)]
    return v

##Trying to write down the rhobar ghost for (non-split) omega^(a-1) + 1 on a-th component
##Probably N=1
def ghost_rhobar_red(p,max_i):
    gg=[0 for i in range(p-1)]
    for a in range(0,p-2,2):
        l1=ab(p,a,p+1-a,0,a,max_i*p)
    	lp=abp(p,a,p+1-a,0,a,max_i*p)
    	g=form_ghost_shell_general(1,p,max_i,a,l1,lp,old=True,no_weight_two=True)
    	gg[a] = g

    return gg

##Trying to write down seqn v with v[k] equal to dim(S_k)_rhobar for 
## rhobar non-split: 0 --> w^beta --> rhobar --> w^alpha --> 0
## 0 <= alpha <= p-2; 1 <= beta <= p-1  (following Serre)
## The formula is: on weights supported by this rhobar the pattern begins with k(rhobar) 1's and then p+1-k(rhobar) 2's. 
## Then add 2 to these p+1 numbers and repeat
## note you need to take M >> 0
def level1_rhobar_dim(p,beta,alpha,M):
    a = min(alpha,beta)
    b = max(alpha,beta)
    krho = 1 + p*a + b
    print "a,b,krho",alpha,beta,krho
    v = [0 for k in range(M)]
    if beta >= alpha:
       c = beta - alpha + 1
    else:
       c = beta - alpha + 1 + p - 1
    for j in range(c):
    	v[krho + j*(p-1)] = 1
    for j in range(c,p+1):
        v[krho + j*(p-1)] = 2
    for j in range(krho+p*(p-1)+1,M):
    	if j%(p-1)==krho%(p-1):
    	   v[j] = v[j-(p+1)*(p-1)]+2
    return v
    
##Trying to write down seqn v with v[k] equal to dim(S_k)_rhobar for 
## rhobar non-split: 0 --> w^beta --> rhobar --> w^alpha --> 0
## 0 <= alpha <= p-2; 1 <= beta <= p-1  (following Serre)
## The formula is: starting at weight beta+alpha+1 2,4,6,8,
def levelp_rhobar_dim(p,beta,alpha,M):
    v = [0 for k in range(M)]
    c = beta + alpha + 1
    for j in range(0,floor(M/(p-1))):
    	if c+j*(p-1)<M:
    	   v[c+j*(p-1)] = 2*(j+1)
    return v
    
##Trying to write down the rhobar ghost for all t-th twists of (non-split) 0 --> w^a --> rhobar --> w^0 --> 0
##Probably N=1 
def ghost_rhobar_red2(p,a,max_i):
    ghosts=[]
    for t in range(0,p-1):
    	print "Trying ",t,a
    	gg=[0 for i in range(p-1)]
    	beta=a+t
	print beta
    	if beta>p-1:
       	   beta-=p-1
    	alpha=t
	print "Calling level 1 with",beta,alpha
    	l1=level1_rhobar_dim(p,beta,alpha,max_i*2*p^2)
    	lp=levelp_rhobar_dim(p,beta,alpha,max_i*2*p^2)
    	## Just guessing here about how many terms to take.
    	newp=[lp[i]-2*l1[i] for i in range(min(len(l1),len(lp)))]
    	g=form_ghost_shell_general(1,p,max_i,(beta+alpha+1)%(p-1),l1,newp,old=True,no_weight_two=True)
    	gg[(beta+alpha+1)%(p-1)] = g
	ghosts+=[gg]

    return ghosts

##single power series attached to omega^a --> rhobar --> 1 non-split
def ghost_rhobar_red_single(p,a,max_i):
    ghosts=[]
    assert a%(p-1)!=1, "That extension not yet programmed"
    beta=a
    if beta>p-1:
       beta-=p-1
    alpha=0
    l1=level1_rhobar_dim(p,beta,alpha,max_i*2*p^2)
    lp=levelp_rhobar_dim(p,beta,alpha,max_i*2*p^2)
    	## Just guessing here about how many terms to take.
    newp=[lp[i]-2*l1[i] for i in range(min(len(l1),len(lp)))]
    g=form_ghost_shell_general(1,p,max_i,(beta+alpha+1)%(p-1),l1,newp,old=True,no_weight_two=True)

    return g

def rim_slopes_rhobar_red(p,max_i):
    g=ghost_rhobar_red(p,max_i)
    v=[wadic_ghost_slopes(g,a,max_i) for a in range(0,p-1,2)]
    print [str([v[a][i+2*p]-v[a][i] for i in range(2*p)]) for a in range(len(v))]
    return [str([v[a][i+1]-v[a][i] for i in range(2*p+10)]) for a in range(len(v))]

def testing(p):
    v=[]
    delta=delta_qexp(30)
    for a in range(0,(p-1)):
        aps=[(q,(delta[q]*q^a)%p) for q in primes(30) if p%q!=0]
	v+=[[counting_fixed_eigensystem(p,1,2,p^2+p+1+40,aps,comp=(12+2*a)%(p-1))]]
    return v

def testing_levelp(p):
    v=[]
    delta=delta_qexp(30)
    for a in range(0,(p-1)):
        aps=[(q,(delta[q]*q^a)%p) for q in primes(30) if p%q!=0]
	v+=[[counting_fixed_levelp_eigensystem(p,1,2,3*p,aps,comp=(12+2*a)%(p-1))]]
    return v

def level1_59_irred(M):
    p=59
    v=[0 for i in range(M)]
    v[16]=1
    for j in range(1,15):
    	v[16+j*(p-1)]=2
    v[16+15*(p-1)]=3
    for j in range(16,p+1):
    	v[16+j*(p-1)]=4
    for j in range(16+p^2-1,M):
        if v[j-(p^2-1)]!=0:
	   v[j]=v[j-(p^2-1)]+4

    return v

def levelp_59_irred(M):
    p=59
    v=[0 for i in range(M)]
    for j in range(floor(M/(p-1))):
    	v[16+j*58]=4*(j+1)
    return v

def ghost_rhobar_59irred(max_i):
    p=59
    l1=level1_59_irred(max_i*2*p^2)
    lp=levelp_59_irred(max_i*2*p^2)
    gg=[0 for i in range(p-1)]
    newp=[lp[i]-2*l1[i] for i in range(min(len(l1),len(lp)))]
    g=form_ghost_shell_general(1,p,max_i,16,l1,newp,old=True,no_weight_two=True)
    gg[16] = g

    return gg