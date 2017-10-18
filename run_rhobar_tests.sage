load("ghost_object.sage")
load("rhobar_tests.sage")

lp=false
verb=false
t=cputime()

print "Testing p=13, N=1, k=12 (Delta)"
print test("DATA/13.1.2.2.3.5.5.7.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=7, N=11, k=2 (X_0(11))"
print test("DATA/7.11.2.5.3.6.5.1.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=11, N=1, k=12 (Delta)"
print test("DATA/11.1.2.9.3.10.5.1.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=23, N=1, k=12 (Delta)"
print test("DATA/23.1.2.22.3.22.5.0.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=7, N=27, k=2 (E=27.a1)"
print test("DATA/7.27.2.0.5.0.11.0.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=59, N=1, k=12 (Delta)"
print test("DATA/59.1.2.35.3.16.5.51.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=19, N=3, k=6 (unique ordinary form)"
print test("DATA/19.3.2.13.5.6.7.17.with_twists.sage",levelp=lp,verbose=verb)

#print

#print "Testing p=13, N=5, k=12 (ordinary form not Delta with mult 2)"
#print test("DATA/13.5.2.2.3.5.7.0.with_twists.sage",levelp=lp,verbose=verb)
## It was Delta


print "Total time",cputime(t)

print "Testing p=23, N=3, k=11 (one of two ordinary forms with irred rhobar and non-trivial nebentype)"
print test("DATA/23.3.2.19.5.10.7.7.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=23, N=3, k=11 (one of two ordinary forms with irred rhobar and non-trivial nebentype)"
print test("DATA/23.3.2.4.5.13.7.7.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=11, N=7, k=3 (irred rhobar and non-trivial nebentype)"
print test("DATA/11.7.2.8.3.0.5.0.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=3, N=11, k=2, X_0(11)[3]"
print test("DATA/3.11.2.1.5.1.7.1.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=3, N=7, k=4, Newform 7.4.1.a"
print test("DATA/3.7.2.2.5.1.11.1.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=3, N=7, k=4, Newform 7.4.1.a"
print test("DATA/3.7.2.2.5.1.11.1.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)

print "Testing p=17, N=2, k=8, E_8"
test_eis("DATA/17.2.3.12.5.11.7.13.with_twists.sage",levelp=lp,verbose=verb)

print "Total time",cputime(t)


