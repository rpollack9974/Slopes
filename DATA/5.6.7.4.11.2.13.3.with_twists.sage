#Rhobar:
#[ [ 7, 4 ], [ 11, 2 ], [ 13, 3 ], [ 17, 4 ], [ 19, 0 ], [ 23, 3 ] ]

#p=5, N=6 arising from weight k=4 (eisenstein)
## argh! it's a twist of E_2!!!

##defect appears to be 9

def dim_tame(k,t):
    p = 5
    t = t % (p-1)
    if k % (p-1) != (4 + 2*t)%(p-1):
        return 0
    if k<4+24:
        if data[k][1][0] == t:
            return len(data[k][1][1])
        elif data[k][2][0] == t:
            return len(data[k][2][1])
        else:
            return "uh oh"
    k0 = (k-4) % 24 + 4
    q = (k-k0)/24
    return 9 * q + dim_tame(k0,t)



data = [ [], [],
[ 2, [ 1,
    []
], [ 3,
    []
] ]
,
[ 3, [ 1,
    []
], [ 3,
    []
] ]
,
[ 4, [ 0,
    [ 0 ]
], [ 2,
    []
] ]
,
[ 5, [ 0,
    []
], [ 2,
    []
] ]
,
[ 6, [ 1,
    []
], [ 3,
    []
] ]
,
[ 7, [ 1,
    []
], [ 3,
    []
] ]
,
[ 8, [ 0,
    [ 0 ]
], [ 2,
    [ 1, 1, 1, 1 ]
] ]
,
[ 9, [ 0,
    []
], [ 2,
    []
] ]
,
[ 10, [ 1,
    [ 1, 1, 1, 1 ]
], [ 3,
    []
] ]
,
[ 11, [ 1,
    []
], [ 3,
    []
] ]
,
[ 12, [ 0,
    [ 0 ]
], [ 2,
    [ 1, 1, 1, 1, 1 ]
] ]
,
[ 13, [ 0,
    []
], [ 2,
    []
] ]
,
[ 14, [ 1,
    [ 1, 1, 1, 1 ]
], [ 3,
    [ 2, 2, 2, 2 ]
] ]
,
[ 15, [ 1,
    []
], [ 3,
    []
] ]
,
[ 16, [ 0,
    [ 0 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1 ]
] ]
,
[ 17, [ 0,
    []
], [ 2,
    []
] ]
,
[ 18, [ 1,
    [ 1, 1, 1, 1 ]
], [ 3,
    [ 2, 2, 2, 2, 2 ]
] ]
,
[ 19, [ 1,
    []
], [ 3,
    []
] ]
,
[ 20, [ 0,
    [ 0, 3, 3, 3, 3 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1 ]
] ]
,
[ 21, [ 0,
    []
], [ 2,
    []
] ]
,
[ 22, [ 1,
    [ 1, 1, 1, 1 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2 ]
] ]
,
[ 23, [ 1,
    []
], [ 3,
    []
] ]
,
[ 24, [ 0,
    [ 0, 3, 3, 3, 3, 3 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1 ]
] ]
,
[ 25, [ 0,
    []
], [ 2,
    []
] ]
,
[ 26, [ 1,
    [ 2, 2, 2, 2, 2, 2, 2, 2 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2 ]
] ]
,
[ 27, [ 1,
    []
], [ 3,
    []
] ]
,
[ 28, [ 0,
    [ 0, 3, 3, 3, 3, 3, 3, 3, 3, 3 ]
], [ 2,
    [ 1, 1, 1, 1, 2, 2, 2, 2, 2 ]
] ]
,
[ 29, [ 0,
    []
], [ 2,
    []
] ]
,
[ 30, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 4 ]
], [ 3,
    [ 3, 3, 3, 3, 3, 3, 3, 3, 3 ]
] ]
,
[ 31, [ 1,
    []
], [ 3,
    []
] ]
,
[ 32, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5 ]
] ]
,
[ 33, [ 0,
    []
], [ 2,
    []
] ]
,
[ 34, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5 ]
], [ 3,
    [ 2, 2, 2, 2, 3, 3, 3, 3, 3 ]
] ]
,
[ 35, [ 1,
    []
], [ 3,
    []
] ]
,
[ 36, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5 ]
] ]
,
[ 37, [ 0,
    []
], [ 2,
    []
] ]
,
[ 38, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 3, 3, 3, 3, 6, 6, 6, 6 ]
] ]
,
[ 39, [ 1,
    []
], [ 3,
    []
] ]
,
[ 40, [ 0,
    [ 0, 3, 3, 3, 3, 4, 4, 4, 4, 4 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6 ]
] ]
,
[ 41, [ 0,
    []
], [ 2,
    []
] ]
,
[ 42, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6 ]
] ]
,
[ 43, [ 1,
    []
], [ 3,
    []
] ]
,
[ 44, [ 0,
    [ 0, 3, 3, 3, 3, 3, 4, 4, 4, 4, 7, 7, 7, 7 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6 ]
] ]
,
[ 45, [ 0,
    []
], [ 2,
    []
] ]
,
[ 46, [ 1,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 6 ]
] ]
,
[ 47, [ 1,
    []
], [ 3,
    []
] ]
,
[ 48, [ 0,
    [ 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 7 ]
], [ 2,
    [ 1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6 ]
] ]
,
[ 49, [ 0,
    []
], [ 2,
    []
] ]
,
[ 50, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8 ]
], [ 3,
    [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5 ]
] ]
,
[ 51, [ 1,
    []
], [ 3,
    []
] ]
,
[ 52, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5 ]
] ]
,
[ 53, [ 0,
    []
], [ 2,
    []
] ]
,
[ 54, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7 ]
], [ 3,
    [ 2, 2, 2, 2, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6 ]
] ]
,
[ 55, [ 1,
    []
], [ 3,
    []
] ]
,
[ 56, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7 ]
] ]
,
[ 57, [ 0,
    []
], [ 2,
    []
] ]
,
[ 58, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 3, 3, 3, 3, 6, 6, 6, 6, 7, 7, 7, 7, 7 ]
] ]
,
[ 59, [ 1,
    []
], [ 3,
    []
] ]
,
[ 60, [ 0,
    [ 0, 3, 3, 3, 3, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 9, 9, 9, 9, 9 ]
] ]
,
[ 61, [ 0,
    []
], [ 2,
    []
] ]
,
[ 62, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 7, 7, 7, 7, 10, 10, 10, 10 ]
] ]
,
[ 63, [ 1,
    []
], [ 3,
    []
] ]
,
[ 64, [ 0,
    [ 0, 3, 3, 3, 3, 3, 4, 4, 4, 4, 7, 7, 7, 7, 8, 8, 8, 8, 8 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 9, 9, 9, 9, 10, 10, 
    10, 10, 10 ]
] ]
,
[ 65, [ 0,
    []
], [ 2,
    []
] ]
,
[ 66, [ 1,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10 ]
] ]
,
[ 67, [ 1,
    []
], [ 3,
    []
] ]
,
[ 68, [ 0,
    [ 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 7, 8, 8, 8, 8, 11, 11, 11, 11 ]
], [ 2,
    [ 1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 9, 9, 9, 9, 10, 10, 
    10, 10, 10 ]
] ]
,
[ 69, [ 0,
    []
], [ 2,
    []
] ]
,
[ 70, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 9, 9, 9, 9, 9 ]
], [ 3,
    [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 
    11, 11, 11, 11 ]
] ]
,
[ 71, [ 1,
    []
], [ 3,
    []
] ]
,
[ 72, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 11, 11, 11, 11, 
    11 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 10, 10, 
    10, 10, 10 ]
] ]
,
[ 73, [ 0,
    []
], [ 2,
    []
] ]
,
[ 74, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 8, 8, 8, 8, 12, 12, 
    12, 12 ]
], [ 3,
    [ 2, 2, 2, 2, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 11, 
    11, 11, 11 ]
] ]
,
[ 75, [ 1,
    []
], [ 3,
    []
] ]
,
[ 76, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 10, 10, 10, 10, 
    10, 10, 10, 10, 10 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9,
    9, 9 ]
] ]
,
[ 77, [ 0,
    []
], [ 2,
    []
] ]
,
[ 78, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 11, 11, 
    11, 11, 11 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 3, 3, 3, 3, 6, 6, 6, 6, 7, 7, 7, 7, 7, 9, 9, 9, 9, 9, 10, 
    10, 10, 10 ]
] ]
,
[ 79, [ 1,
    []
], [ 3,
    []
] ]
,
[ 80, [ 0,
    [ 0, 3, 3, 3, 3, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 
    10, 10, 10, 10, 10 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 12, 12, 12, 12 ]
] ]
,
[ 81, [ 0,
    []
], [ 2,
    []
] ]
,
[ 82, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 11, 11, 
    11, 11, 11, 11, 11, 11, 11 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 7, 7, 7, 7, 10, 10, 10, 10, 10, 
    10, 10, 10, 10 ]
] ]
,
[ 83, [ 1,
    []
], [ 3,
    []
] ]
,
[ 84, [ 0,
    [ 0, 3, 3, 3, 3, 3, 4, 4, 4, 4, 7, 7, 7, 7, 8, 8, 8, 8, 8, 11, 11, 11, 11, 
    11, 11, 11, 11, 11 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 9, 9, 9, 9, 10, 10, 
    10, 10, 10, 12, 12, 12, 12, 12 ]
] ]
,
[ 85, [ 0,
    []
], [ 2,
    []
] ]
,
[ 86, [ 1,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 12, 12, 
    12, 12, 12, 12, 12, 12, 12 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 
    12, 12, 12, 12, 12, 12, 12, 12 ]
] ]
,
[ 87, [ 1,
    []
], [ 3,
    []
] ]
,
[ 88, [ 0,
    [ 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 7, 8, 8, 8, 8, 11, 11, 11, 11, 
    12, 12, 12, 12, 12 ]
], [ 2,
    [ 1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 9, 9, 9, 9, 10, 10, 
    10, 10, 10, 13, 13, 13, 13, 13, 13, 13, 13, 13 ]
] ]
,
[ 89, [ 0,
    []
], [ 2,
    []
] ]
,
[ 90, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 9, 9, 9, 9, 9, 13, 13, 
    13, 13, 13, 13, 13, 13, 13 ]
], [ 3,
    [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 
    11, 11, 11, 11, 14, 14, 14, 14, 14 ]
] ]
,
[ 91, [ 1,
    []
], [ 3,
    []
] ]
,
[ 92, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 11, 11, 11, 11, 
    11, 12, 12, 12, 12, 15, 15, 15, 15 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 10, 10, 
    10, 10, 10, 14, 14, 14, 14, 14, 14, 14, 14, 14 ]
] ]
,
[ 93, [ 0,
    []
], [ 2,
    []
] ]
,
[ 94, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 8, 8, 8, 8, 12, 12, 
    12, 12, 13, 13, 13, 13, 13 ]
], [ 3,
    [ 2, 2, 2, 2, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 11, 
    11, 11, 11, 14, 14, 14, 14, 15, 15, 15, 15, 15 ]
] ]
,
[ 95, [ 1,
    []
], [ 3,
    []
] ]
,
[ 96, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 10, 10, 10, 10, 
    10, 10, 10, 10, 10, 15, 15, 15, 15, 15 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9,
    9, 9, 14, 14, 14, 14, 14, 14, 14, 14, 14 ]
] ]
,
[ 97, [ 0,
    []
], [ 2,
    []
] ]
,
[ 98, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 11, 11, 
    11, 11, 11, 12, 12, 12, 12, 16, 16, 16, 16 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 3, 3, 3, 3, 6, 6, 6, 6, 7, 7, 7, 7, 7, 9, 9, 9, 9, 9, 10, 
    10, 10, 10, 13, 13, 13, 13, 15, 15, 15, 15, 15 ]
] ]
,
[ 99, [ 1,
    []
], [ 3,
    []
] ]
,
[ 100, [ 0,
    [ 0, 3, 3, 3, 3, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 
    10, 10, 10, 10, 10, 14, 14, 14, 14, 14, 16, 16, 16, 16 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 12, 12, 12, 12, 13, 13, 13, 13, 13 ]
] ]
,
[ 101, [ 0,
    []
], [ 2,
    []
] ]
,
[ 102, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 11, 11, 
    11, 11, 11, 11, 11, 11, 11, 15, 15, 15, 15, 15 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 7, 7, 7, 7, 10, 10, 10, 10, 10, 
    10, 10, 10, 10, 13, 13, 13, 13, 14, 14, 14, 14, 14 ]
] ]
,
[ 103, [ 1,
    []
], [ 3,
    []
] ]
,
[ 104, [ 0,
    [ 0, 3, 3, 3, 3, 3, 4, 4, 4, 4, 7, 7, 7, 7, 8, 8, 8, 8, 8, 11, 11, 11, 11, 
    11, 11, 11, 11, 11, 14, 14, 14, 14, 14, 15, 15, 15, 15 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 9, 9, 9, 9, 10, 10, 
    10, 10, 10, 12, 12, 12, 12, 12, 13, 13, 13, 13, 16, 16, 16, 16 ]
] ]
,
[ 105, [ 0,
    []
], [ 2,
    []
] ]
,
[ 106, [ 1,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 12, 12, 
    12, 12, 12, 12, 12, 12, 12, 15, 15, 15, 15, 15, 15, 15, 15, 15 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 
    12, 12, 12, 12, 12, 12, 12, 12, 14, 14, 14, 14, 14 ]
] ]
,
[ 107, [ 1,
    []
], [ 3,
    []
] ]
,
[ 108, [ 0,
    [ 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 7, 8, 8, 8, 8, 11, 11, 11, 11, 
    12, 12, 12, 12, 12, 14, 14, 14, 14, 14, 15, 15, 15, 15 ]
], [ 2,
    [ 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 9, 9, 9, 9, 10, 10, 
    10, 10, 10, 13, 13, 13, 13, 13, 13, 13, 13, 13, 16, 16, 16, 16, 16 ]
] ]
,
[ 109, [ 0,
    []
], [ 2,
    []
] ]
,
[ 110, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 9, 9, 9, 9, 9, 13, 13, 
    13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 15, 15, 15, 15 ]
], [ 3,
    [ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 10, 10, 10, 10, 
    11, 11, 11, 11, 14, 14, 14, 14, 14, 14, 14, 14, 14, 17, 17, 17, 17 ]
] ]
,
[ 111, [ 1,
    []
], [ 3,
    []
] ]
,
[ 112, [ 0,
    [ 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 11, 11, 11, 11, 
    11, 12, 12, 12, 12, 15, 15, 15, 15, 15, 15, 15, 15, 15 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 7, 7, 10, 10, 
    10, 10, 10, 14, 14, 14, 14, 14, 14, 14, 14, 14, 16, 16, 16, 16, 16, 16, 16, 
    16, 16 ]
] ]
,
[ 113, [ 0,
    []
], [ 2,
    []
] ]
,
[ 114, [ 1,
    [ 1, 1, 1, 1, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 12, 12, 
    12, 12, 13, 13, 13, 13, 13, 16, 16, 16, 16, 16, 16, 16, 16, 16 ]
], [ 3,
    [ 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 11, 
    11, 11, 11, 14, 14, 14, 14, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17 ]
] ]
,
[ 115, [ 1,
    []
], [ 3,
    []
] ]
,
[ 116, [ 0,
    [ 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 17, 17, 17 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8,
    8, 8, 14, 14, 14, 14, 14, 14, 14, 14, 14, 17, 17, 17, 17, 17, 17, 17, 17, 17
    ]
] ]
,
[ 117, [ 0,
    []
], [ 2,
    []
] ]
,
[ 118, [ 1,
    [ 1, 1, 1, 1, 5, 5, 5, 5, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 
    10, 10, 10, 11, 11, 11, 11, 16, 16, 16, 16, 17, 17, 17, 17, 17 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 4, 4, 4, 4, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9,
    9, 9, 12, 12, 12, 12, 15, 15, 15, 15, 15, 18, 18, 18, 18, 18, 18, 18, 18, 18
    ]
] ]
,
[ 119, [ 1,
    []
], [ 3,
    []
] ]
,
[ 120, [ 0,
    [ 0, 3, 3, 3, 3, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 13, 13, 13, 13, 13, 16, 16, 16, 16, 19, 19, 19, 19, 19 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 6, 6, 6, 6, 7, 7, 7, 7, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 11, 11, 11, 11, 12, 12, 12, 12, 12, 18, 18, 18, 18, 18, 18, 18, 18, 18
    ]
] ]
,
[ 121, [ 0,
    []
], [ 2,
    []
] ]
,
[ 122, [ 1,
    [ 1, 1, 1, 1, 5, 5, 5, 5, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 14, 14, 14, 14, 14, 15, 15, 15, 15, 20, 
    20, 20, 20 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 7, 7, 7, 7, 7, 8, 8, 8, 8, 10, 10, 10, 10, 10, 
    10, 10, 10, 10, 12, 12, 12, 12, 13, 13, 13, 13, 13, 19, 19, 19, 19, 19, 19, 
    19, 19, 19 ]
] ]
,
[ 123, [ 1,
    []
], [ 3,
    []
] ]
,
[ 124, [ 0,
    [ 0, 3, 3, 3, 3, 3, 5, 5, 5, 5, 8, 8, 8, 8, 9, 9, 9, 9, 9, 11, 11, 11, 11, 
    11, 11, 11, 11, 11, 13, 13, 13, 13, 13, 14, 14, 14, 14, 17, 17, 17, 17, 20, 
    20, 20, 20, 20 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 6, 6, 6, 6, 7, 7, 7, 7, 10, 10, 10, 10, 11, 
    11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 15, 15, 15, 15, 16, 16, 
    16, 16, 16 ]
] ]
,
[ 125, [ 0,
    []
], [ 2,
    []
] ]
,
[ 126, [ 1,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 14, 14, 14, 14, 14, 14, 14, 14, 14, 18, 
    18, 18, 18, 18 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 7, 7, 7, 7, 7, 7, 7, 7, 7, 11, 11, 11, 11, 11, 
    12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 17, 17, 17, 17, 17, 17, 
    17, 17, 17 ]
] ]
,
[ 127, [ 1,
    []
], [ 3,
    []
] ]
,
[ 128, [ 0,
    [ 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 8, 8, 8, 8, 8, 9, 9, 9, 9, 12, 12, 12, 12, 
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 17, 17, 17, 17, 18, 
    18, 18, 18, 18 ]
], [ 2,
    [ 1, 1, 1, 1, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 7, 7, 7, 7, 10, 10, 10, 10, 11, 
    11, 11, 11, 11, 13, 13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 16, 
    16, 16, 16, 19, 19, 19, 19 ]
] ]
,
[ 129, [ 0,
    []
], [ 2,
    []
] ]
,
[ 130, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 4, 5, 5, 5, 5, 9, 9, 9, 9, 10, 10, 10, 10, 10, 14,
    14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 18, 18, 
    18, 18, 18, 19, 19, 19, 19 ]
], [ 3,
    [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 11, 11, 11, 11, 11, 
    12, 12, 12, 12, 14, 14, 14, 14, 14, 14, 14, 14, 14, 16, 16, 16, 16, 17, 17, 
    17, 17, 17 ]
] ]
,
[ 131, [ 1,
    []
], [ 3,
    []
] ]
,
[ 132, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 12, 12, 12, 12, 
    12, 13, 13, 13, 13, 15, 15, 15, 15, 15, 15, 15, 15, 15, 17, 17, 17, 17, 18, 
    18, 18, 18, 18 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, 6, 6, 6, 6, 6, 9, 9, 9, 9, 11, 11, 
    11, 11, 11, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 
    15, 15, 19, 19, 19, 19, 19 ]
] ]
,
[ 133, [ 0,
    []
], [ 2,
    []
] ]
,
[ 134, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 9, 9, 9, 9, 13, 13, 
    13, 13, 14, 14, 14, 14, 14, 16, 16, 16, 16, 16, 16, 16, 16, 16, 18, 18, 18, 
    18, 18, 19, 19, 19, 19 ]
], [ 3,
    [ 2, 2, 2, 2, 3, 3, 3, 3, 3, 7, 7, 7, 7, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10, 
    12, 12, 12, 12, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 
    17, 17, 17, 20, 20, 20, 20 ]
] ]
,
[ 135, [ 1,
    []
], [ 3,
    []
] ]
,
[ 136, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 11, 11, 11, 11, 
    11, 11, 11, 11, 11, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 18, 
    18, 18, 18, 18 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 
    10, 10, 10, 15, 15, 15, 15, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 17, 17, 
    17, 17, 19, 19, 19, 19, 19, 19, 19, 19, 19 ]
] ]
,
[ 137, [ 0,
    []
], [ 2,
    []
] ]
,
[ 138, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 12, 12, 
    12, 12, 12, 13, 13, 13, 13, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 
    18, 18, 19, 19, 19, 19 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 3, 3, 3, 3, 6, 6, 6, 6, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 
    11, 11, 11, 11, 14, 14, 14, 14, 16, 16, 16, 16, 16, 18, 18, 18, 18, 18, 18, 
    18, 18, 18, 20, 20, 20, 20, 20 ]
] ]
,
[ 139, [ 1,
    []
], [ 3,
    []
] ]
,
[ 140, [ 0,
    [ 0, 3, 3, 3, 3, 4, 4, 4, 4, 4, 9, 9, 9, 9, 9, 9, 9, 9, 9, 11, 11, 11, 11, 
    11, 11, 11, 11, 11, 15, 15, 15, 15, 15, 17, 17, 17, 17, 19, 19, 19, 19, 19, 
    19, 19, 19, 19, 21, 21, 21, 21 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 10, 10, 10, 10, 10, 
    10, 10, 10, 10, 13, 13, 13, 13, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 
    19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19 ]
] ]
,
[ 141, [ 0,
    []
], [ 2,
    []
] ]
,
[ 142, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 16, 16, 16, 16, 16, 17, 17, 17, 17, 20, 
    20, 20, 20, 20, 20, 20, 20, 20 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 8, 8, 8, 8, 11, 11, 11, 11, 11, 
    11, 11, 11, 11, 14, 14, 14, 14, 15, 15, 15, 15, 15, 20, 20, 20, 20, 20, 20, 
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20 ]
] ]
,
[ 143, [ 1,
    []
], [ 3,
    []
] ]
,
[ 144, [ 0,
    [ 0, 3, 3, 3, 3, 3, 4, 4, 4, 4, 7, 7, 7, 7, 9, 9, 9, 9, 9, 12, 12, 12, 12, 
    12, 12, 12, 12, 12, 15, 15, 15, 15, 15, 16, 16, 16, 16, 19, 19, 19, 19, 21, 
    21, 21, 21, 21, 21, 21, 21, 21, 21 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 10, 10, 10, 10, 11, 
    11, 11, 11, 11, 13, 13, 13, 13, 13, 14, 14, 14, 14, 17, 17, 17, 17, 18, 18, 
    18, 18, 18, 21, 21, 21, 21, 21, 21, 21, 21, 21 ]
] ]
,
[ 145, [ 0,
    []
], [ 2,
    []
] ]
,
[ 146, [ 1,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 16, 16, 16, 16, 16, 16, 16, 16, 16, 20, 
    20, 20, 20, 20, 22, 22, 22, 22, 22, 22, 22, 22 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 6, 11, 11, 11, 11, 11, 
    13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 19, 19, 19, 19, 19, 19, 
    19, 19, 19, 22, 22, 22, 22, 22, 22, 22, 22, 22 ]
] ]
,
[ 147, [ 1,
    []
], [ 3,
    []
] ]
,
[ 148, [ 0,
    [ 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 7, 9, 9, 9, 9, 12, 12, 12, 12, 
    13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 16, 16, 16, 16, 19, 19, 19, 19, 20, 
    20, 20, 20, 20, 23, 23, 23, 23, 23, 23, 23, 23, 23 ]
], [ 2,
    [ 1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 10, 10, 10, 10, 11, 
    11, 11, 11, 11, 14, 14, 14, 14, 14, 14, 14, 14, 14, 17, 17, 17, 17, 17, 18, 
    18, 18, 18, 21, 21, 21, 21, 22, 22, 22, 22, 22 ]
] ]
,
[ 149, [ 0,
    []
], [ 2,
    []
] ]
,
[ 150, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 10, 10, 10, 10, 10, 14,
    14, 14, 14, 14, 14, 14, 14, 14, 16, 16, 16, 16, 16, 16, 16, 16, 16, 20, 20, 
    20, 20, 20, 21, 21, 21, 21, 24, 24, 24, 24, 24 ]
], [ 3,
    [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 11, 11, 11, 11, 11, 
    12, 12, 12, 12, 15, 15, 15, 15, 15, 15, 15, 15, 15, 18, 18, 18, 18, 19, 19, 
    19, 19, 19, 23, 23, 23, 23, 23, 23, 23, 23, 23 ]
] ]
,
[ 151, [ 1,
    []
], [ 3,
    []
] ]
,
[ 152, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 12, 12, 12, 12, 
    12, 13, 13, 13, 13, 16, 16, 16, 16, 16, 16, 16, 16, 16, 19, 19, 19, 19, 20, 
    20, 20, 20, 20, 24, 24, 24, 24, 24, 24, 24, 24, 24 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 9, 9, 9, 9, 11, 11, 
    11, 11, 11, 15, 15, 15, 15, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 17, 17, 
    17, 17, 21, 21, 21, 21, 21, 22, 22, 22, 22, 25, 25, 25, 25 ]
] ]
,
[ 153, [ 0,
    []
], [ 2,
    []
] ]
,
[ 154, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 9, 9, 9, 9, 13, 13, 
    13, 13, 14, 14, 14, 14, 14, 17, 17, 17, 17, 17, 17, 17, 17, 17, 20, 20, 20, 
    20, 20, 21, 21, 21, 21, 24, 24, 24, 24, 25, 25, 25, 25, 25 ]
], [ 3,
    [ 2, 2, 2, 2, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 
    12, 12, 12, 12, 15, 15, 15, 15, 16, 16, 16, 16, 16, 18, 18, 18, 18, 18, 19, 
    19, 19, 19, 22, 22, 22, 22, 23, 23, 23, 23, 23 ]
] ]
,
[ 155, [ 1,
    []
], [ 3,
    []
] ]
,
[ 156, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 11, 11, 11, 11, 
    11, 11, 11, 11, 11, 16, 16, 16, 16, 16, 18, 18, 18, 18, 18, 18, 18, 18, 20, 
    20, 20, 20, 20, 24, 24, 24, 24, 24, 24, 24, 24, 24 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 10, 10, 
    10, 10, 10, 15, 15, 15, 15, 15, 15, 15, 15, 15, 18, 18, 18, 18, 18, 18, 18, 
    18, 18, 21, 21, 21, 21, 21, 21, 21, 21, 21, 25, 25, 25, 25, 25 ]
] ]
,
[ 157, [ 0,
    []
], [ 2,
    []
] ]
,
[ 158, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 12, 12, 
    12, 12, 12, 13, 13, 13, 13, 17, 17, 17, 17, 18, 18, 18, 18, 18, 20, 20, 20, 
    20, 20, 21, 21, 21, 21, 24, 24, 24, 24, 25, 25, 25, 25, 25 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 3, 3, 3, 3, 6, 6, 6, 6, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10, 
    11, 11, 11, 11, 14, 14, 14, 14, 16, 16, 16, 16, 16, 19, 19, 19, 19, 19, 19, 
    19, 19, 19, 22, 22, 22, 22, 22, 23, 23, 23, 23, 26, 26, 26, 26 ]
] ]
,
[ 159, [ 1,
    []
], [ 3,
    []
] ]
,
[ 160, [ 0,
    [ 0, 3, 3, 3, 3, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 11, 11, 11, 11, 
    11, 11, 11, 11, 11, 15, 15, 15, 15, 15, 17, 17, 17, 17, 20, 20, 20, 20, 20, 
    20, 20, 20, 20, 23, 23, 23, 23, 24, 24, 24, 24, 24 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 9, 9, 9, 9, 9, 10, 
    10, 10, 10, 13, 13, 13, 13, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 
    19, 19, 21, 21, 21, 21, 21, 21, 21, 21, 21, 25, 25, 25, 25, 25, 26, 26, 26, 
    26 ]
] ]
,
[ 161, [ 0,
    []
], [ 2,
    []
] ]
,
[ 162, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 12, 12, 
    12, 12, 12, 12, 12, 12, 12, 16, 16, 16, 16, 16, 17, 17, 17, 17, 21, 21, 21, 
    21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 25, 25, 25, 25, 25 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 7, 7, 7, 7, 10, 10, 10, 10, 11, 
    11, 11, 11, 11, 14, 14, 14, 14, 15, 15, 15, 15, 15, 20, 20, 20, 20, 20, 20, 
    20, 20, 20, 22, 22, 22, 22, 22, 22, 22, 22, 22, 26, 26, 26, 26, 26 ]
] ]
,
[ 163, [ 1,
    []
], [ 3,
    []
] ]
,
[ 164, [ 0,
    [ 0, 3, 3, 3, 3, 3, 4, 4, 4, 4, 7, 7, 7, 7, 8, 8, 8, 8, 8, 12, 12, 12, 12, 
    12, 12, 12, 12, 12, 15, 15, 15, 15, 15, 16, 16, 16, 16, 19, 19, 19, 19, 21, 
    21, 21, 21, 21, 23, 23, 23, 23, 23, 24, 24, 24, 24, 27, 27, 27, 27 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 6, 6, 6, 6, 9, 9, 9, 9, 10, 10, 
    10, 10, 10, 13, 13, 13, 13, 13, 14, 14, 14, 14, 17, 17, 17, 17, 18, 18, 18, 
    18, 18, 22, 22, 22, 22, 22, 22, 22, 22, 22, 25, 25, 25, 25, 25, 26, 26, 26, 
    26 ]
] ]
,
[ 165, [ 0,
    []
], [ 2,
    []
] ]
,
[ 166, [ 1,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 13, 13, 
    13, 13, 13, 13, 13, 13, 13, 16, 16, 16, 16, 16, 16, 16, 16, 16, 20, 20, 20, 
    20, 20, 23, 23, 23, 23, 23, 23, 23, 23, 25, 25, 25, 25, 25 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 
    13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 19, 19, 19, 19, 19, 19, 
    19, 19, 19, 23, 23, 23, 23, 23, 23, 23, 23, 23, 26, 26, 26, 26, 26, 26, 26, 
    26, 26 ]
] ]
,
[ 167, [ 1,
    []
], [ 3,
    []
] ]
,
[ 168, [ 0,
    [ 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 7, 8, 8, 8, 8, 11, 11, 11, 11, 
    13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 16, 16, 16, 16, 19, 19, 19, 19, 20, 
    20, 20, 20, 20, 24, 24, 24, 24, 24, 24, 24, 24, 24, 27, 27, 27, 27, 27 ]
], [ 2,
    [ 1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 9, 9, 9, 9, 10, 10, 
    10, 10, 10, 14, 14, 14, 14, 14, 14, 14, 14, 14, 17, 17, 17, 17, 17, 18, 18, 
    18, 18, 21, 21, 21, 21, 22, 22, 22, 22, 22, 25, 25, 25, 25, 25, 26, 26, 26, 
    26 ]
] ]
,
[ 169, [ 0,
    []
], [ 2,
    []
] ]
,
[ 170, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 9, 9, 9, 9, 9, 14, 14, 
    14, 14, 14, 14, 14, 14, 14, 16, 16, 16, 16, 16, 16, 16, 16, 16, 20, 20, 20, 
    20, 20, 21, 21, 21, 21, 24, 24, 24, 24, 25, 25, 25, 25, 25, 28, 28, 28, 28 ]
], [ 3,
    [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 
    11, 11, 11, 11, 15, 15, 15, 15, 15, 15, 15, 15, 15, 18, 18, 18, 18, 19, 19, 
    19, 19, 19, 23, 23, 23, 23, 23, 23, 23, 23, 23, 26, 26, 26, 26, 26, 26, 26, 
    26, 26 ]
] ]
,
[ 171, [ 1,
    []
], [ 3,
    []
] ]
,
[ 172, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 11, 11, 11, 11, 
    11, 13, 13, 13, 13, 16, 16, 16, 16, 16, 16, 16, 16, 16, 19, 19, 19, 19, 20, 
    20, 20, 20, 20, 24, 24, 24, 24, 24, 24, 24, 24, 24, 27, 27, 27, 27, 27, 27, 
    27, 27, 27 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 10, 10, 
    10, 10, 10, 15, 15, 15, 15, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 17, 17, 
    17, 17, 21, 21, 21, 21, 21, 22, 22, 22, 22, 25, 25, 25, 25, 25, 26, 26, 26, 
    26 ]
] ]
,
[ 173, [ 0,
    []
], [ 2,
    []
] ]
,
[ 174, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 8, 8, 8, 8, 12, 12, 
    12, 12, 14, 14, 14, 14, 14, 17, 17, 17, 17, 17, 17, 17, 17, 17, 20, 20, 20, 
    20, 20, 21, 21, 21, 21, 24, 24, 24, 24, 25, 25, 25, 25, 25, 28, 28, 28, 28, 
    28 ]
], [ 3,
    [ 2, 2, 2, 2, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 11, 
    11, 11, 11, 15, 15, 15, 15, 16, 16, 16, 16, 16, 18, 18, 18, 18, 18, 19, 19, 
    19, 19, 22, 22, 22, 22, 23, 23, 23, 23, 23, 26, 26, 26, 26, 26, 26, 26, 26, 
    26 ]
] ]
,
[ 175, [ 1,
    []
], [ 3,
    []
] ]
,
[ 176, [ 0,
    [ 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 10, 10, 10, 10, 
    10, 10, 10, 10, 10, 16, 16, 16, 16, 16, 18, 18, 18, 18, 18, 18, 18, 18, 20, 
    20, 20, 20, 20, 24, 24, 24, 24, 24, 24, 24, 24, 24, 27, 27, 27, 27, 27, 27, 
    27, 27, 27 ]
], [ 2,
    [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9,
    9, 9, 15, 15, 15, 15, 15, 15, 15, 15, 15, 18, 18, 18, 18, 18, 18, 18, 18, 
    18, 21, 21, 21, 21, 21, 21, 21, 21, 21, 25, 25, 25, 25, 25, 27, 27, 27, 27, 
    27, 27, 27, 27 ]
] ]
,
[ 177, [ 0,
    []
], [ 2,
    []
] ]
,
[ 178, [ 1,
    [ 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 11, 11, 
    11, 11, 11, 13, 13, 13, 13, 17, 17, 17, 17, 18, 18, 18, 18, 18, 20, 20, 20, 
    20, 20, 21, 21, 21, 21, 24, 24, 24, 24, 25, 25, 25, 25, 25, 28, 28, 28, 28, 
    28, 28, 28, 28, 28 ]
], [ 3,
    [ 2, 2, 2, 2, 2, 3, 3, 3, 3, 6, 6, 6, 6, 7, 7, 7, 7, 7, 9, 9, 9, 9, 9, 10, 
    10, 10, 10, 14, 14, 14, 14, 16, 16, 16, 16, 16, 19, 19, 19, 19, 19, 19, 19, 
    19, 19, 22, 22, 22, 22, 22, 23, 23, 23, 23, 26, 26, 26, 26, 27, 27, 27, 27, 
    27 ]
] ]
,
[ 179, [ 1,
    []
], [ 3,
    []
] ]
]