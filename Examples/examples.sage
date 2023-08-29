load("../phts_hyp.sage")

# Example for local heights code.
# Tested using Sage 9.8
# The final examples take very long.


print("We test a few odd degree examples in genus 2 and compare our code to Balakrishnan-Besser.")
p = 13
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f = x*(x-1)*(x-2)*(x-5)*(x-6)
C = HyperellipticCurve(f)
P = C(3,6)
Q = C(3,-6)
R = C(10,120)
S = C(10,-120)
T = C(0,0)
t0 = cputime()
hpqrs = height_four_affine_points(P, Q, R, S, prec)
t1 = cputime()
hrspq = height_four_affine_points(R, S, P, Q, prec)
t2 = cputime()
hpqst = height_four_affine_points(P, Q, S, T, prec)
t3 = cputime()
hpqrt = height_four_affine_points(P, Q, R, T, prec)
t4 = cputime()
print("individual heights", hpqrs, hrspq)
print("check symmetry", hpqrs-hrspq)
print("check path independence",hpqrs+hpqst-hpqrt)
print("times for heights using our algorithm", t1-t0, t2-t1, t3-t2, t4-t3)
oldhpqrs = C.height([ (1, P), (-1, Q)], [ (1, R), (-1, S)])
t5 = cputime()
print("check equality with Balakrishnan-Besser", hpqrs-oldhpqrs)
print("time for height using Balakrishnan-Besser", t5-t4)



print("\nNow, one genus 3 example.")
p = 11
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f = x^7 + (x-1)^2*(x-4)^2*(x-9)*(x-16)
C = HyperellipticCurve(f)
P = C(1,1)
Q = C(4,2^7)
R = C(9,3^7)
S = C(16,4^7)
T = C(0,48)
t0 = cputime()
hpqrs = height_four_affine_points(P, Q, R, S, prec)
t1 = cputime()
hrspq = height_four_affine_points(R, S, P, Q, prec)
t2 = cputime()
hpqst = height_four_affine_points(P, Q, S, T, prec)
t3 = cputime()
hpqrt = height_four_affine_points(P, Q, R, T, prec)
t4 = cputime()
oldhpqrs = C.height([ (1, P), (-1, Q)], [ (1, R), (-1, S)])
print("check path independence", hpqrs+hpqst-hpqrt)
print("times for heights using our algorithm", t1-t0, t2-t1, t3-t2, t4-t3)
t5 = cputime()
oldhpqrt = C.height([ (1, P), (-1, Q)], [ (1, R), (-1, T)])
t6 = cputime()
print("check equality between our algorithm and Balakrishnan-Besser",hpqrt-oldhpqrt)
print("time for height using Balakrishnan-Besser", t5-t4)

print("\nIn this even degree example, of genus 2, we check the model independence of h_p, by testing h_p(\\infty_- - \\infty_+, P-Q) on C, and on the curves X obtained by sending x to 1/x, using both our algorithm and Balakrishnan-Besser. We also show that precomputing psi(x^gdx/y) can save a lot of time.")
p = 5
prec = 10
K = pAdicField(p, prec)
x = polygen(K)
C = HyperellipticCurve(x^6+2*x^2+x)
P = C(1,2)
NP = opposite_affine_point(P)
Q = C.lift_x(6)
NQ = opposite_affine_point(Q)
R = C.lift_x(2)
#NR = opposite_affine_point(R)
#S = C.lift_x(7)
#NS = opposite_affine_point(S)
#T = C.lift_x(2016)
#NT = opposite_affine_point(T)
X = HyperellipticCurve(x^5+2*x^4+1)
A = X(0,-1)
NA = X(0,1)
B = X(1/(P[0]),(P[1])/(P[0]^3))
NB = opposite_affine_point(B)
D = X(1/(Q[0]),(Q[1])/(Q[0]^3))
ND = opposite_affine_point(D)
#E = X(1/(R[0]),(R[1])/(R[0]^3))
#NE = opposite_affine_point(E)
#F = X(1/(S[0]),(S[1])/(S[0]^3))
#NF = opposite_affine_point(F)
#G = X(1/(T[0]),(T[1])/(T[0]^3))
#NG = opposite_affine_point(G)
t0 = cputime()
hPC = height_infinities(P, NP, prec)
t1 = cputime()
M = psi_omega_for_infinities(C, prec)
t2 = cputime()
hPC2 = height_infinities(P, NP, prec, M=M)
t3 = cputime()
hPX = height_four_affine_points(A, NA, B, NB, prec)
print("height same on C and X?", hPC - hPX)
hPXBB = X.height([ (1, A), (-1, NA)], [ (1, B), (-1, NB)])
print("height on X agrees with BB12?", hPX - hPXBB)
print("There seems to a small precision loss in BB12 in some cases.")
hQC = height_infinities(Q, NQ, prec)
hQX = height_four_affine_points(A, NA, D, ND, prec)
print("height same on C and X?", hPC - hPX)
hQXBB = X.height([ (1, A), (-1, NA)], [ (1, D), (-1, ND)])
print("height on X agrees with BB12?", hPX - hPXBB)
print("There seems to a small precision loss in BB12 in some cases.")
#height_infinities(R, NR, prec)
#height_four_affine_points(A, NA, E, NE, prec)
#X.height([ (1, A), (-1, NA)], [ (1, E), (-1, NE)])
#height_infinities(S, NS, prec)
#height_four_affine_points(A, NA, F, NF, prec)
#X.height([ (1, A), (-1, NA)], [ (1, F), (-1, NF)])
#height_infinities(T, NT, prec)
#height_four_affine_points(A, NA, G, NG, prec)
#X.height([ (1, A), (-1, NA)], [ (1, G), (-1, NG)])
print("time for height without precomputed psi(x^gdx/2y)", t1-t0)
print("time to compute psi(x^gdx/2y)", t2-t1)
print("time for height with precomputed psi(x^gdx/2y)", t3-t2)

print("\n Now similarly, in genus 3.")
p = 11
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f1 = x^8+3*x^5+4*x^4+x
C = HyperellipticCurve(f1)
P = C(1,3)
Q = C(-1,1)
f2 = x^7+4*x^4+3*x^3+1
X = HyperellipticCurve(f2)
R = X(0,-1)
S = X(0,1)
P2 = X(1,3)
Q2 = X(-1,1)
t0 = cputime()
hC = height_infinities(P, Q, prec)
t1 = cputime()
hX = height_four_affine_points(R, S, P2, Q2, prec)
t2 = cputime()
hXBB = X.height([ (1, R), (-1, S)], [ (1, P2), (-1, Q2)])
t3 = cputime()
print("time to compute height on C using diffs of infinite pts", t1-t0)
print("time to compute height on X using 4 affine pts", t2-t1)
print("time to compute height on X using BB12", t3-t2)



print("\n We show that our algorithm can deal with different complementary subspaces W; namely the unit root subspace and the one obtained from a symplectic basis. We see that both heights are symmetric.")
p = 11
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f = (x^4+1)^2+x*(x-1)*(x-2)*(x-3)
C = HyperellipticCurve(f)
P = C(0,1)
Q = C(1,2)
R = C(2, 17)
S = C(3, 82)
t0 = cputime()
h0 = height_four_affine_points(P, Q, R, S, prec)
print("height wrt default choice of subspace", h0)
t1 = cputime()
t1-t0
W = changebasistounitroot(C)
h1 = height_four_affine_points(P, Q, R, S, prec, BasisW = W)
print("height wrt unit of subspace", h1)
print("As expected, since unit root subspace is default choice")
W2 = changebasistosymplectic(C)
h2 = height_four_affine_points(P, Q, R, S, prec, BasisW = W2)
print("height wrt symplectic subspace", h2)
h3 = height_four_affine_points(R, S, P, Q, prec, BasisW = W2)
print("height wrt symplectic subspace symmetric?", h2-h3)


print("\n We consider heights of the type hp(P-w(P),Q-w(Q)), where we allow one of them to be in a residue disc of a point at infinity. The resulting height is symmetric")
p = 13
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f = x*(x-1)*(x-2)*(x-3)*(x-4)*(x-5)+1
C = HyperellipticCurve(f)
P = C(0,1)
Q = C(1,-1)
R = C(2,1)
S = C(3,-1)
T = C(4,1)
U = C(5,-1)
SamplePts = [Q,R,S,T,U]
ListPtsNegVal = [C.lift_x(i+1/p) for i in range(5)]
for i in range(5):
    h1 = height_both_antisymmetric(SamplePts[i], ListPtsNegVal[i], prec)
    h2 = height_both_antisymmetric(ListPtsNegVal[i], SamplePts[i], prec)
    print("symmetric?", h1-h2)

print("\n We consider heights on both odd and even degree hyperelliptic curves, and allow points to be in bad discs, i.e., Weiestrass and discs at infinity, in various combinations. Symmetry still holds. In the odd degree examples, we get the same result as BB12.")
p = 11
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f1 = x*(x*(x-1)*(x-4)*(x-9)*(x-16)+1)
C = HyperellipticCurve(f1)
P = C(1,1)
Q = C(4,2)
R = C.lift_x(1/p)
S = C(9,3)
T = C(16,4)
f2 = x^5+(1-x)*(1-4*x)*(1-9*x)*(1-16*x)
X = HyperellipticCurve(f2)
P2 = X(1/P[0],P[1]/P[0]^3)
Q2 = X(1/Q[0],Q[1]/Q[0]^3)
R2 = X(1/R[0],R[1]/R[0]^3)
S2 = X(1/S[0],S[1]/S[0]^3)
T2 = X(1/T[0],T[1]/T[0]^3)
t0 = cputime()
h1 = height_four_affine_points(P, Q, R, S, prec)
t1 = cputime()
h2 = height_four_affine_points(R, S, P, Q, prec)
print("points in Weierstrass discs and infinite discs in even degree: heights symmetric?", h1-h2)
t2 = cputime()
h3 = height_four_affine_points(P2, Q2, R2, S2, prec)
t3 = cputime()
h4 = height_four_affine_points(R2, S2, P2, Q2, prec)
print("points in Weierstrass discs and infinite disc in odd degree: heights symmetric?", h3-h4)
t4 = cputime()
h5 = height_four_affine_points(P, Q, S, T, prec)
t5 = cputime()
h6 = height_four_affine_points(S, T, P, Q, prec)
t6 = cputime()
print("points in Weierstrass discs in even degree: heights symmetric?", h5-h6)
h7 = height_four_affine_points(P2, Q2, S2, T2, prec)
t7 = cputime()
h8 = height_four_affine_points(S2, T2, P2, Q2, prec)
print("points in Weierstrass discs in odd degree: heights symmetric?", h7-h8)
t8 = cputime()
h9 = height_four_affine_points(P, Q, R, T, prec)
t9 = cputime()
h10 = height_four_affine_points(R, T, P, Q, prec)
t10 = cputime()
print("points in Weierstrass discs and infinite discs in even degree: heights symmetric?", h9-h10)
h11 = height_four_affine_points(P2, Q2, R2, T2, prec)
t11 = cputime()
h12 = height_four_affine_points(R2, T2, P2, Q2, prec)
t12 = cputime()
print("points in Weierstrass discs and infinite disc in odd degree: heights symmetric?", h11-h12)
h3BB = X.height([ (1, P2), (-1, Q2)], [ (1, R2), (-1, S2)])
t13 = cputime()
h7BB = X.height([ (1, P2), (-1, Q2)], [ (1, S2), (-1, T2)])
t14 = cputime()
h11BB = X.height([ (1, P2), (-1, Q2)], [ (1, R2), (-1, T2)])
t15 = cputime()
print("points in Weierstrass discs and infinite disc in odd degree: same as BB12?", h3-h3BB)
print("time for our algorithm", t3-t2)
print("time for BB12", t13-t12)
print("points in Weierstrass discs in odd degree: same as BB12?", h7-h7BB)
print("time for our algorithm", t7-t6)
print("time for BB12", t14-t13)
print("points in Weierstrass discs and infinite disc in odd degree: same as BB12?", h11-h11BB)
print("time for our algorithm", t11-t10)
print("time for BB12", t15-t14)


print("\n We test the function hp(\\infty_- -Q,R-S) for two genus 2 curves related by x<->1/x.")
p = 11
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f = x*(x-1)*(x-2)*(x-3)*(x-4)*(x-5)+1
C = HyperellipticCurve(f)
P = C(0,1)
L1 = [C(i,1) for i in range(1,6)]
f2 = (1-x)*(1-2*x)*(1-3*x)*(1-4*x)*(1-5*x)+x^6
X = HyperellipticCurve(f2)
Q = X(0,-1)
L2 = [X(1/P[0],P[1]/P[0]^3) for P in L1]
h1 = [height_infinity_minus_three_affine_points(L1[i], L1[j], L1[k], prec) for i in range(5) for j in range(i+1,5) for k in range(j+1,5)]
h2 = [height_four_affine_points(Q, L2[i], L2[j], L2[k], prec) for i in range(5) for j in range(i+1,5) for k in range(j+1,5)]
print("heights model independent?", [h1[i]-h2[i] for i in range(len(h1))])


print("\n We test both functions hp(\\infty_- -Q,R-S) and hp(\\infty_+ -Q,R-S) for two genus 3 curves related by x<->1/x.")
p = 11
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f = x^8-8*x^7+6*x^6-3*x^5+7*x^4+2*x^3-x^2-10*x+1
C = HyperellipticCurve(f)
P = C(0,1)
Q = C.lift_x(-45)
R = C.lift_x(3)
S = C.lift_x(5)
L1 = [Q, R, S]
f2 = x^8-10*x^7-x^6+2*x^5+7*x^4-3*x^3+6*x^2-8*x+1
X = HyperellipticCurve(f2)
P2 = X(0,-1)
L2 = [X(1/P[0],P[1]/P[0]^4) for P in L1]
print("height model independent?", height_infinity_minus_three_affine_points(Q, P, R, prec) - height_infinity_plus_three_affine_points(L2[1], P2, L2[0], prec))
print("height model independent?", height_infinity_minus_three_affine_points(Q, P, S, prec) - height_infinity_plus_three_affine_points(L2[2], P2, L2[0], prec))
print("height model independent?", height_infinity_minus_three_affine_points(R, P, S, prec) - height_infinity_plus_three_affine_points(L2[2], P2, L2[1], prec))
print("height model independent?", height_four_affine_points(R, S, P, Q, prec) - height_infinity_plus_three_affine_points(L2[0], L2[1], L2[2], prec))
print("height model and path independent?", height_infinity_minus_three_affine_points(Q, P, R, prec) - height_infinity_minus_three_affine_points(Q, P, S, prec) - height_four_affine_points(P2, L2[0], L2[2], L2[1], prec))


print("\n We check our specialized code for 'tiny heights' hp(\\infty_- -\\infty_+,P-Q) for P and Q in the same disc against the generic code.")
p = 11
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f = x^6-3*x^5+7*x^4+2*x^3-x^2-10*x+1
C = HyperellipticCurve(f)
P = C.lift_x(0)
Q = C.lift_x(6)
R = C.lift_x(9)
P2 = C.lift_x(11)
Q2 = C.lift_x(-5)
R2 = C.lift_x(-2)
t0 = cputime()
print("height computed using tiny code", tiny_height_infinities(P, P2, prec))
print("height computed using generic code", height_infinities(P, P2, prec))
print("height computed using tiny code", tiny_height_infinities(Q, Q2, prec))
print("height computed using generic code", height_infinities(Q, Q2, prec))
print("height computed using tiny code", tiny_height_infinities(R, R2, prec))
print("height computed using generic code", height_infinities(R, R2, prec))



print("\n We test the function height_first_antisymmetric(P, R, S, prec), which computes hp(P-w(P),R-S).")
p = 11
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f = x^6+(p*x-1)*(x^4-1)
C = HyperellipticCurve(f)
P = C(1/p,-1/p^3)
NP = C(1/p,1/p^3)
R = C(1,1)
NR = C(1,-1)
print("The following four heights should all agree, as they all should be equal to hp(R-NR, P-NP)")
height_both_antisymmetric_infinity_disc(R, P, prec)
height_first_antisymmetric(P, R, NR, prec)
height_first_antisymmetric(R, P, NP, prec)
height_four_affine_points(P, NP, R, NR, prec)


print("\n WARNING: The following examples take very long to check!")



print("\n Now we test the function height_divisors(div1, div2, prec) where we plug in two degree zero divisors satisfying the conditions of our algorithm and we see that the order of decompositions of divisors does not matter, as expected.")
p = 7
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
f = x^6-x^5+x^2+1
C = HyperellipticCurve(f)
P1 = C.lift_x(0)
Q1 = C.lift_x(1)
R1 = C.lift_x(2)
S1 = C.lift_x(4)
P2 = C.lift_x(7)
Q2 = C.lift_x(8)
R2 = C.lift_x(9)
S2 = C.lift_x(11)
P3 = C.lift_x(14)
Q3 = C.lift_x(15)
P = [P1, P2, P3]
Q = [Q1, Q2, Q3]
R = [R1, R2]
S = [S1, S2]
D1 = pair_out(P, Q)
D2 = pair_out(R, S)
hd1d2 = height_divisors(D1, D2, prec)
print("hp(D1, D2), computed using height_divisors:", hd1d2)
print("Now compute hp(D1, D2) using different choices of orders of decompositions")
permutation_indices1 = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
permutation_indices2 = [[0, 1], [1, 0]]
M = psi_omega_for_infinities(C, prec)
for I1 in permutation_indices1:
    for I2 in permutation_indices1:
        for I3 in permutation_indices2:
            for I4 in permutation_indices2:
                print("height equal to hp(D1, D2)?", sum(height_four_affine_points(P[I1[i]], Q[I2[i]], R[I3[j]], S[I4[j]], prec) for i in range(3) for j in range(2)) - hd1d2)



print("\n In this multipurpose example, we test more cases where points are in bad discs. Furthermore, we test the function hp(\\infty-Q,R-S) on odd degree models. Results agree, as expected.")
p = 13
prec = 7
K = pAdicField(p, prec)
x = polygen(K)
g = 3
f1 = x*(x*(x-1)*(x-4)*(x-9)*(x-16)*(x-25)*(x-36)+1)
C = HyperellipticCurve(f1)
P = C(0,0)
L1 = [C(i^2,i) for i in range(1,5)]
L2 = [C.lift_x(i^2*p^2) for i in range(1,3)]
L3 = [C.lift_x(i+1/p) for i in range(2)]
L4 = L1 + L2 + L3
f2 = x^7+(1-x)*(1-4*x)*(1-9*x)*(1-16*x)*(1-25*x)*(1-36*x)
X = HyperellipticCurve(f2)
K1 = [X(1/L1[i][0],L1[i][1]/L1[i][0]^(g+1)) for i in range(4)]
K2 = [X(1/L2[i][0],L2[i][1]/L2[i][0]^(g+1)) for i in range(2)]
K3 = [X(1/L3[i][0],L3[i][1]/L3[i][0]^(g+1)) for i in range(2)]
K4 = K1 + K2 + K3
M = psi_omega_for_infinities(C, prec)
heights1 = [height_four_affine_points(L1[i], L1[j], L4[k], L4[l], prec) for i in range(4) for j in range(i+1,4) for k in range(j+1,8) for l in range(k+1,8)]
numberofexamples = len(heights1)
heights2 = [height_four_affine_points(L4[k], L4[l], L1[i], L1[j], prec) for i in range(4) for j in range(i+1,4) for k in range(j+1,8) for l in range(k+1,8)]
print("heights symmetric?", [heights1[i] - heights2[i] for i in range(numberofexamples)])
heights3 = [height_four_affine_points(K1[i], K1[j], K4[k], K4[l], prec) for i in range(4) for j in range(i+1,4) for k in range(j+1,8) for l in range(k+1,8)]
heights4 = [height_four_affine_points(K4[k], K4[l], K1[i], K1[j], prec) for i in range(4) for j in range(i+1,4) for k in range(j+1,8) for l in range(k+1,8)]
print("heights symmetric?", [heights3[i] - heights4[i] for i in range(numberofexamples)])
print("heights model-independent?",[heights1[i] - heights3[i] for i in range(numberofexamples)])
heights5 = [height_four_affine_points(P, L4[k], L1[i], L1[j], prec) for i in range(4) for j in range(i+1,4) for k in range(j+1,8)]
heights6 = [height_four_affine_points(L1[i], L1[j], P, L4[k], prec) for i in range(4) for j in range(i+1,4) for k in range(j+1,8)]
numberofexamples2 = len(heights5)
print("heights symmetric?", [heights5[i] - heights6[i] for i in range(numberofexamples2)])
heights7 = [height_infinity_three_affine_points(K4[k], K1[i], K1[j], prec) for i in range(4) for j in range(i+1,4) for k in range(j+1,8)]
print("heights model independent?", [heights5[i] - heights7[i] for i in range(numberofexamples2)])
