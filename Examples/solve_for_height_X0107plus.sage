
# Solve for the global p-adic height as a bilinar pairing on X0+(107)
# for p=7; used for quadratic Chabauty + Mordell-Weil sieve
# to compute the rational points on X0+(107)

t0 = cputime()
load("X0107plus_qc_setup.m")
load("../phts_hyp.sage")
load("X0107plus_p7_init.sage")
h11 = height_divisors(pair_out(D11, E11), pair_out(F11, G11), 10)
h12 = height_divisors(pair_out(D12, E12), pair_out(F12, G12), 10)
h22 = height_divisors(pair_out(D22, E22), pair_out(F22, G22), 10)
print("h11", h11)
print("h12", h12)
print("h22", h22)
t1 = cputime()
print("Time", t1-t0)

