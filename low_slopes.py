def write_low_slopes_to_file(p,N,num_coefs,min_weight,max_weight):
	g = ghost(p,N,num_coefs=num_coefs)
	filename = "DATA/low_slopes." + str(p) + "." + str(N) + "." + str(num_coefs) 
	for k in range(min_weight,max_weight,2):
		f = open(filename,"a")
		slopes = g.slopes(k)
		f.write(str(k) + ":" + str(slopes) + "\n")
		f.close()
	return "done"
