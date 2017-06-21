
##data comes from MAGMA program finding_modp_repns.magma except the data needs to massaged by
##removing all *s and inserting two entries in the start of the list so that data[k] is a list of 
##length 3: the first entry is k, the second entry is of the form [t,[list of slopes]] where t is the 
##smaller twist so that the rhobar slopes occur in weight k.  The third entry has the form
##[t +(p-1/2),[list of slopes]] giving the slopes which occur in the other possible twist
##
##rb comes from rbdata class (defined in recursive_formulas.sage)
##terms is the number of coefficients the ghost series is computed withw
def test_rhobar_data(data,rb,terms=10,verbose=false):
	p = rb.p
	g = [ghost(rbdata=rb,twist=t,terms=terms) for t in range(0,p-1)]
	for k in range(2,len(data)):
		if k%2 == 0:
			t = min(((k-rb.krbar)/2) % (p-1),((k-rb.krbar)/2 + (p-1)/2) % (p-1))
			if verbose:
				print "Weight",k,"and twist",t
				print "True slopes:", data[k][1][1]
				print "Ghost slopes:", g[t].slopes(k)[0:S(k,t,rb)]
				print "Weight",k,"and twist",(t+(p-1)/2) % (p-1)
				print "True slopes:", data[k][2][1]
				print "Ghost slopes:", g[(t+(p-1)/2)%(p-1)].slopes(k)[0:S(k,(t+(p-1)/2)%(p-1),rb)]
			if data[k][1][1] != g[t].slopes(k)[0:S(k,t,rb)]:
				print "Error in weight",k,"and twist",t
				print "True slopes:", data[k][1][1]
				print "Ghost slopes:", g[t].slopes(k)[0:S(k,t,rb)]
			if data[k][2][1] != g[(t+(p-1)/2)%(p-1)].slopes(k)[0:S(k,(t+(p-1)/2)%(p-1),rb)]:
				print "Error in weight",k,"and twist",(t+(p-1)/2) % (p-1)
				print "True slopes:", data[k][2][1]
				print "Ghost slopes:", g[(t+(p-1)/2)%(p-1)].slopes(k)[0:S(k,(t+(p-1)/2)%(p-1),rb)]
