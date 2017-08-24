print "Testing p=13, N=1, k=12 (Delta)"
print test("DATA/13.1.2.2.3.5.5.7.with_twists.sage",levelp=false)

print

print "Testing p=7, N=11, k=2 (X_0(11))"
print test("DATA/7.11.2.5.3.6.5.1.with_twists.sage",levelp=false)

print

print "Testing p=11, N=1, k=12 (Delta)"
print test("DATA/11.1.2.9.3.10.5.1.with_twists.sage",levelp=false)

print

print "Testing p=23, N=1, k=12 (Delta)"
print test("DATA/23.1.2.22.3.22.5.0.with_twists.sage",levelp=false)

print

print "Testing p=7, N=7, k=2 (E=27.a1)"
print test("DATA/7.27.2.0.5.0.11.0.with_twists.sage",levelp=false)
