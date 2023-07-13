"""
Check the p-adic Birch and Swinnerton-Dyer conjecture from
Balakrishnan-Muller-Stein for J0(67)+ and the good ordinary primes 
p = 11, 29, 31, 71, 89.

This is a long computation, mostly due to the computation of the p-adic
L-values, which is done using overconvergent modular symbols. Change prec
below to make it faster.
"""

from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
from sage.modular.pollack_stevens.padic_lseries import pAdicLseries

def Lseries(space, p, prec = 5):
    """
    Compute p-adic L-series of a modular symbol space corresponding to
    a  modular abelian surface and the corresponding p-adic multipliers.
    """
    phi = ps_modsym_from_simple_modsym_space(space)
    ap = phi.Tq_eigenvalue(p,prec)
    c1,c2 = phi.completions(p,prec)
    phi1,psi1 = c1
    phi2,psi2 = c2
    phi1p = phi1.p_stabilize_and_lift(p,ap = psi1(ap), M = prec) 
    L1 = pAdicLseries(phi1p)                                     
    eps1 = L1.interpolation_factor(psi1(ap))
    phi2p = phi2.p_stabilize_and_lift(p,ap = psi2(ap), M = prec) 
    L2  = pAdicLseries(phi2p)  
    eps2 = L1.interpolation_factor(psi2(ap))
    return [L1,L2], [eps1,eps2]

def Lvaluerank2(space, p, prec = 5):
    """
    Compute special value of p-adic L-series of modular symbol space of
    weight 2, corresponding to a modular abelian surface of
    rank 2 and p-adic multiplier.
    """
    Ls, epss = Lseries(space, p, prec=prec)
    assert Ls[0][0] == 0  
    assert Ls[1][0] == 0  
    return Ls[0][1]*Ls[1][1], epss[0]*epss[1]

spaces = ModularSymbols(67,2,1).cuspidal_submodule().new_subspace().decomposition()
space = spaces[1] # This corresponds to X0+(67) (e.g. by comparing with BMS16)

# Find good representatives for regulator computations
load("X067plus_regulator_setup.m") 
delta_plus = 1/4 # see Table 4.4 in BMS16

print("\n-----------\n")

p = 11
print("Check BSD for p = ", p)
try:
    prec = 12
    Lval1, eps1 = Lvaluerank2(space, p, prec = prec) 
    load("X067plus_p11_init.sage")
    load("../phts_hyp.sage")
    h11 = height_divisors(pair_out(D11, E11), pair_out(F11, G11), prec)-away11
    h22 = height_divisors(pair_out(D22, E22), pair_out(F22, G22), prec)-away22
    h12 = height_divisors(pair_out(D12, E12), pair_out(F12, G12), prec)-away12
    reg = h11*h22/(m11*n11*m22*n22) - h12^2/(m12*n12)^2
    print("11-adic regulator", reg)
    K = Qp(p, prec)
    # All integral BSD quantities are 1, see Table 4.2 in BMS16
    bsd_quot = Lval1*delta_plus/eps1 / (reg/log(K(1+p))^2)
    print("BSD quotient", bsd_quot)
except NotImplementedError:
    print("Error: p needs to split in RM field.")

print("\n-----------\n")

p = 29
print("Check BSD for p = ", p)
try:
    prec = 11
    Lval1, eps1 = Lvaluerank2(space, p, prec = prec) 
    load("X067plus_p29_init.sage")
    load("../phts_hyp.sage")
    h11 = height_divisors(pair_out(D11, E11), pair_out(F11, G11), prec)-away11
    h22 = height_divisors(pair_out(D22, E22), pair_out(F22, G22), prec)-away22
    h12 = height_divisors(pair_out(D12, E12), pair_out(F12, G12), prec)-away12
    reg = h11*h22/(m11*n11*m22*n22) - h12^2/(m12*n12)^2
    print("29-adic regulator", reg)
    K = Qp(p, prec)
    # All integral BSD quantities are 1, see Table 4.2 in BMS16
    bsd_quot = Lval1*delta_plus/eps1 / (reg/log(K(1+p))^2)
    print("BSD quotient", bsd_quot)
except NotImplementedError:
    print("Error: p needs to split in RM field.")

print("\n-----------\n")

p = 31
print("Check BSD for p = ", p)
try:
    prec = 11
    Lval1, eps1 = Lvaluerank2(space, p, prec = prec) 
    load("X067plus_p31_init.sage")
    h11 = height_divisors(pair_out(D11, E11), pair_out(F11, G11), prec)-away11
    h22 = height_divisors(pair_out(D22, E22), pair_out(F22, G22), prec)-away22
    h12 = height_divisors(pair_out(D12, E12), pair_out(F12, G12), prec)-away12
    reg = h11*h22/(m11*n11*m22*n22) - h12^2/(m12*n12)^2
    print("31-adic regulator", reg)
    K = Qp(p, prec)
    # All integral BSD quantities are 1, see Table 4.2 in BMS16
    bsd_quot = Lval1*delta_plus/eps1 / (reg/log(K(1+p))^2)
    print("BSD quotient", bsd_quot)
except NotImplementedError:
    print("Error: p needs to split in RM field.")

print("\n-----------\n")


p = 71
print("Check BSD for p = ", p)
try:
    prec = 12
    Lval1, eps1 = Lvaluerank2(space, p, prec = prec) 
    load("X067plus_p71_init.sage")
    h11 = height_divisors(pair_out(D11, E11), pair_out(F11, G11), prec)-away11
    h22 = height_divisors(pair_out(D22, E22), pair_out(F22, G22), prec)-away22
    h12 = height_divisors(pair_out(D12, E12), pair_out(F12, G12), prec)-away12
    reg = h11*h22/(m11*n11*m22*n22) - h12^2/(m12*n12)^2
    print("71-adic regulator", reg)
    K = Qp(p, prec)
    # All integral BSD quantities are 1, see Table 4.2 in BMS16
    bsd_quot = Lval1*delta_plus/eps1 / (reg/log(K(1+p))^2)
    print("BSD quotient", bsd_quot)
except NotImplementedError:
    print("Error: p needs to split in RM field.")

print("\n-----------\n")

p = 89
print("Check BSD for p = ", p)
try:
    prec = 11
    Lval1, eps1 = Lvaluerank2(space, p, prec = prec) 
    load("X067plus_p89_init.sage")
    h11 = height_divisors(pair_out(D11, E11), pair_out(F11, G11), prec)-away11
    h22 = height_divisors(pair_out(D22, E22), pair_out(F22, G22), prec)-away22
    h12 = height_divisors(pair_out(D12, E12), pair_out(F12, G12), prec)-away12
    reg = h11*h22/(m11*n11*m22*n22) - h12^2/(m12*n12)^2
    print("89-adic regulator", reg)
    K = Qp(p, prec)
    # All integral BSD quantities are 1, see Table 4.2 in BMS16
    bsd_quot = Lval1*delta_plus/eps1 / (reg/log(K(1+p))^2)
    print("BSD quotient", bsd_quot)
except NotImplementedError:
    print("Error: p needs to split in RM field.")

print("\n-----------\n")
