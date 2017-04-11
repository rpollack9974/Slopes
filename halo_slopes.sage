from sage.geometry.newton_polygon import NewtonPolygon

def slopes_on_rim(p):
    G = DirichletGroup(p^2)
    psi = G.0
    chi = psi^(p-1)
    K = psi.base_ring()
    OK = K.ring_of_integers()
    print "Factoring ",p    
    pp = K.prime_above(p)
    e = pp.ramification_index()
    omega = psi^p
    print "Finding omega"
    r = 0
    vals = []
    mv = 0
    while mv == 0:
        vals = []    
        r = r+1
        for a in range(1,p):
            vals += [(omega(a)^r-a).valuation(pp)]
        mv = min(vals)

    omega = omega^r
	
    seeds= []
    for j in range(0,(p-1)/2):
        print "Computing the ",j,"-th modular symbol space"
        M = ModularSymbols(chi*omega^(2*j),2,1,base_ring=K).cuspidal_subspace()
        print "---computing Hecke"
        f = M.hecke_polynomial(p)
        f = f.coefficients()
        f.reverse()
        v = []
        for r in range(len(f)):
            v += [(r,f[r].valuation(pp)/e)]
        ans = NewtonPolygon(v).slopes()
        print j,ans
        seeds += [ans]

    all_slopes = []
    for comp in range(0,p-1):
        if comp % 2 == 0:
            sl = []
            for a in range(1,(p+1)/2):
                sl += [a*(p-1)] #eisenstein stuff
            for j in range(0,(p-1)/2):
                i = (comp-2-2*j) % (p-1)
                i = ZZ(i/2)
                for s in range(0,len(seeds[i])):
                    sl += [(seeds[i][s]+j)*(p-1)]
                sl.sort()
            all_slopes += [[comp,sl]]

    return all_slopes
