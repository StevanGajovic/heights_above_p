
load("../phts_hyp.sage")


p = 7
prec = 10
K = pAdicField(p, prec)
t = polygen(K)
C = HyperellipticCurve(t^5+(t-1)*(t-4)*(t-9)*(t-16))
P = C(1,1)
Q = C(4,2^5)
R = C(9,3^5)
S = C(16,4^5)
t1 = cputime()
h1 = height_four_affine_points(P, Q, R, S, 10)
t2 = cputime()
h2  = C.height([ (1, P), (-1, Q)], [ (1, R), (-1, S)])
t3 = cputime()
print("Difference between heights")
print(h1-h2)
print("time for our algorithm:", t2-t1)
#time for our algorithm: 2.315206
print("time for Balakrishnan-Besser:", t3-t2)
#time for Balakrishnan-Besser: 7.239461


#genus 3
p = 11
prec = 10
K = pAdicField(p, prec)
x = polygen(K)
f = x^7 + (x-1)^2*(x-4)^2*(x-9)*(x-16)
C = HyperellipticCurve(f)
P = C(1,1)
Q = C(4,2^7)
R = C(9,3^7)
S = C(16,4^7)
T = C(0,48)
t1 = cputime()
h1 = height_four_affine_points(P, Q, R, S, prec)
t2 = cputime()
h2  = C.height([ (1, P), (-1, Q)], [ (1, R), (-1, S)])
t3 = cputime()
print("Difference between heights")
print(h1-h2)
print("time for our algorithm:", t2-t1)
#time for our algorithm: 6.312056999999999
print("time for Balakrishnan-Besser:", t3-t2)
#time for Balakrishnan-Besser: 28.003038999999998



# genus 4
p = 23
prec = 20
K = pAdicField(p, prec)
t = polygen(K)
C = HyperellipticCurve(t^9+(t-1)^2*(t-4)^2*(t-9)^2*(t-16)^2)
P = C(1,1)
Q = C(4,2^9)
R = C(9,3^9)
S = C(16,4^9)
t1 = cputime()
h1 = height_four_affine_points(P, Q, R, S, 20)
t2 = cputime()
h2 = C.height([ (1, P), (-1, Q)], [ (1, R), (-1, S)])
t3 = cputime()
print("Difference between heights")
print(h1-h2)
print("time for our algorithm:", t2-t1)
#time for our algorithm: 113.74753900000002
print("time for Balakrishnan-Besser:", t3-t2)
#time for Balakrishnan-Besser: 2749.718002



# large prime, genus 2
p = 503
prec = 10
K = pAdicField(p, prec)
t = polygen(K)
C = HyperellipticCurve(t^5+(t-1)*(t-4)*(t-9)*(t-16))
P = C(1,1)
Q = C(4,2^5)
R = C(9,3^5)
S = C(16,4^5)
t1 = cputime()
h3 = height_four_affine_points(P, Q, R, S, prec)
t2 = cputime()
print("time for our algorithm:", t2-t1)
#time for our algorithm: 246.00850600000012


# large precision, genus 2
p = 7
prec = 300
K = pAdicField(p, prec)
t = polygen(K)
C = HyperellipticCurve(t^5+(t-1)*(t-4)*(t-9)*(t-16))
P = C(1,1)
Q = C(4,2^5)
R = C(9,3^5)
S = C(16,4^5)
t1 = cputime()
h2 = height_four_affine_points(P, Q, R, S, prec)
t2 = cputime()
print("time for our algorithm:", t2-t1)
#time for our algorithm: 648.1203089999999




#
# genus 17
p = 11
prec = 7
K = pAdicField(p, prec)
t = polygen(K)
C = HyperellipticCurve(t^35+(t-1)^2*(t-4)^2*(t-9)^2*(t-16)^2)
P = C(1,1)
Q = C(4,2^35)
R = C(9,3^35)
S = C(16,4^35)
t1 = cputime()
h1 = height_four_affine_points(P, Q, R, S, prec)
t2 = cputime()
print("time for our algorithm:", t2-t1)
#time for our algorithm: 4299.250149
#genus 2
