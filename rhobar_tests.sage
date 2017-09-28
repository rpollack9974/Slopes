load("ghost_object.sage")
load("recursive_dimension_formulas.sage")

def level_N_to_level_Np(slopes,k,t,rb,f1s=None,f2s=None):
	ans = slopes[:]
	assert len(slopes)==0 or max(slopes) <= (k-2)/2,"Slopes too high"
	if f1s == None:
		dp = Sp(k,t,rb)
	else:
		dp = f2s[t](k) + 2*f1s[t](k)
	d = len(slopes)
	dnew = dp - 2*d
	for j in range(dnew):
		ans += [(k-2)/2]
	for j in range(d-1,-1,-1):
		ans += [k-1-slopes[j]]
	return ans

def enhance_data_to_levelp(data2,rb,f1s=None,f2s=None):
#	data = deepcopy(data2)
	for m in range(2,len(data)):
		t=data[m][1][0]
		data[m][1][1] = level_N_to_level_Np(data[m][1][1],m,t,rb,f1s=f1s,f2s=f2s)
		t=data[m][2][0]
		data[m][2][1] = level_N_to_level_Np(data[m][2][1],m,t,rb,f1s=f1s,f2s=f2s)



##data comes from MAGMA program finding_modp_repns.magma except the data needs to massaged by
##removing all *s and inserting two entries in the start of the list so that the k-th entry
##corresponds to weight k.  Then data[k] is a list of 
##length 3: the first entry is k, the second entry is of the form [t,[list of slopes]] where t is the 
##smaller twist so that the rhobar slopes occur in weight k.  The third entry has the form
##[t +(p-1/2),[list of slopes]] giving the slopes which occur in the other possible twist
##
##rb comes from rbdata class (defined in recursive_formulas.sage)
##terms is the number of coefficients the ghost series is computed with
def test_rhobar_data(data,rb,f1s=None,f2s=None,verbose=false,levelp=False):
	p = rb.p
	A = max([len(data[k][1][1]) for k in range(2,len(data))])
	B = max([len(data[k][2][1]) for k in range(2,len(data))])
	terms = 2*max(A,B)+1
	if verbose:
		print "Requires",terms,"terms!"
	g = []
	for t in range(0,p-1):
		if verbose:
			print "forming ghost with twist",t
		if f1s == None:
			g += [ghost(rbdata=rb,twist=t,terms=terms)]
		else:
			g += [ghost(terms=terms,p=p,f1=f1s[t],f2=f2s[t],comp=(rb.krbar+2*t)%(p-1))]
	passed = True
	for k in range(2,len(data)):
		if len(data[k][1][1]) != 0:
			t = min(((k-rb.krbar)/2) % (p-1),((k-rb.krbar)/2 + (p-1)/2) % (p-1))
			tp = (t + (p-1)/2)%(p-1)
			data1 = data[k][1][1]
			data2 = data[k][2][1]
			if levelp == True:
				if f1s == None:
					gdata1 = g[t].slopes(k,terms=Sp(k,t,rb))
				else:
					print "using",f2s[t](k) + 2*f1s[t](k),"terms"
					gdata1 = g[t].slopes(k,terms=f2s[t](k) + 2*f1s[t](k))
			else:
				if f1s == None:
					gdata1 = g[t].slopes(k,terms=S(k,t,rb))
				else:
					gdata1 = g[t].slopes(k,terms=f1s[t](k))

			if levelp == True:
				if f1s == None:
					gdata2 = g[tp].slopes(k,terms=Sp(k,tp,rb))
				else:
					gdata2 = g[tp].slopes(k,terms=f2s[tp](k) + 2*f1s[tp](k))
			else:
				if f1s == None:
					gdata2 = g[tp].slopes(k,terms=S(k,tp,rb))
				else:
					gdata2 = g[tp].slopes(k,terms=f1s[tp](k))

			if verbose:
				print "Weight",k,"and twist",t
				print "True slopes:", data1
				print "Ghost slopes:", gdata1
				print "Weight",k,"and twist",tp
				print "True slopes:", data2
				print "Ghost slopes:", gdata2
			if data1 != gdata1:
				passed = False
				print "*********************************************************************************"
				print "Error in weight",k,"and twist",t
				print "True slopes:", data1
				print "Ghost slopes:", gdata1
			if data2 != gdata2:
				passed = False
				print "*********************************************************************************"
				print "Error in weight",k,"and twist",tp
				print "True slopes:", data2
				print "Ghost slopes:", gdata2
	return passed


def test(filename,verbose=false,levelp=false):
	load(filename)
	print "Rhobar given by",rb
	if levelp:
		enhance_data_to_levelp(data,rb)     
		print "Working at level Np up to weight",len(data)-1
	else:
		print "Working at level N up to weight",len(data)-1
	return test_rhobar_data(data,rb,verbose=verbose,levelp=levelp) 	

def test_eis(filename,verbose=false,levelp=false):
	load(filename)
	print "Rhobar given by",rb
	if levelp:
		enhance_data_to_levelp(data,rb)     
		print "Working at level Np up to weight",len(data)-1
	else:
		print "Working at level N up to weight",len(data)-1
	return test_rhobar_data(data,rb,verbose=verbose,levelp=levelp,f1s=f1s,f2s=f2s) 	