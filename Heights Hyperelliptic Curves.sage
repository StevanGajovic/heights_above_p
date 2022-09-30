"""
This is a code to compute the local compontent of the p-adic (Coleman-Gross) height above p on hyperelliptic curve $C:y^2=f(x)$, i.e., for given points $P,Q,R,S\in C(\mathbb{Q}_p)$ satisfying conditions below, computes $h_p(P-Q,R-S)$:
Conditions on points $P,Q,R,S\in C(\mathbb{Q}_p)$:
(1) Either all four points are affine, or one pair of $\{P,Q\}$ and $\{R,S\}$ is the difference $\infty_- -\infty_+$ of points at infinity.
(2) Residue discs of $P$ and $Q$ are distinct of residue discs of $R$, $\iota(R)$, $S$, $\iota(S)$, where $\iota$ is the hyperelliptic involution.
(3) Affine points do not belong to residue discs at infinity, unless $\{P,Q\}$ or $\{R,S\}$ are in the same residue disc.
Then:
height_infinities(C, P, Q, prec):    Computes h_p(\infty_- - \infty_+, P-Q) up to precision prec
height_four_affine_points(C, P, Q, R, S, prec):    Computes h_p(P-Q, R-S) up to precision prec
"""

def local_coordinates_at_infinity_even(H, prec = 20, name = 't'):
    """
    For an even degree model, only gives local coordinate above one point
    """
    g = H.genus()
    pol1,pol2 = H.hyperelliptic_polynomials()
    K = LaurentSeriesRing(H.base_ring(), name)
    t = K.gen()
    K.set_default_prec(prec+2)
    L = PolynomialRing(K,'x')
    x = L.gen()
    i = 0
    x = t**-1
    f = pol1(x)
    y = (pol1.leading_coefficient()).sqrt()*t**-(g+1)
    for i in range((RR(log(prec)/log(2))).ceil()):
        y = (y + f/y)/2
    return (x+O(t**(prec+2)) , y+O(t**(prec+2)))

def h1dr_basis(self):
    """
    Extract a basis of H^1_dR(C) from H^1_MW(U)
    """
    g = self.genus()
    f = self.hyperelliptic_polynomials()[0]
    if f.degree() == 2*g+1:
        x,y = self.monsky_washnitzer_gens()
        w = self.invariant_differential()
        return [x^i*w for i in range(2*g)]    
    if f.degree() == 2*g+2:
        x,y = local_coordinates_at_infinity_even(self)
        wg = (x^g)*x.derivative()/(2*y)
        c = wg.residue()
        L = []
        for i in range(g+1,2*g+1):
            diff = x^i*x.derivative()/(2*y)
            r = diff.residue()
            xx,yy = self.monsky_washnitzer_gens()
            w = self.invariant_differential()
            wg_xxyy = xx^(g)*w
            diff_xxyy = xx^i*w
            L = L + [(diff_xxyy-r/c *wg_xxyy)]
        return [xx^i*w for i in range(g)] + L

def matrix_of_frobenius_on_h1dr(C, p, prec):
    """
    Used only for an even degree model, computes the matrix of Fronbenius action on H^1_dR(C) from Frobenius action on H^1_MW(U) (doing some linear algebra)
    """
    g = C.genus()
    basis = h1dr_basis(C)
    X = C.change_ring(QQ)
    F = X.matrix_of_frobenius(p,prec)
    FH = matrix(C.base_ring(), 2*g,2*g)
    for i in range(g):
        for j in range(g):
            FH[i,j] = F[i,j]
    for i in range(g):
        for j in range(g,2*g):
            FH[i,j] = F[i,j+1] + basis[j].coeffs()[0][0][g]*F[i,g]
    for i in range(g,2*g):
        for j in range(g):
            FH[i,j] = F[i+1,j]
    for i in range(g,2*g):
        for j in range(g,2*g):
            FH[i,j] = F[i+1,j+1] + basis[j].coeffs()[0][0][g]*F[i+1,g]
    return FH

def cup_product_matrix_even_basis(self, basis):
    """
    For an even degree model, computes the cup product matrix on H^1_dR(C)
    """
    g = self.genus()
    #Warning: Sage sees only one point above infinity
    x,y = local_coordinates_at_infinity_even(self)
    xprime = x.derivative()
    w = [0 for i in range(2*g)]
    for i in range(2*g):
        coeffs = basis[i].coeffs()[0][0]  # x-coefficients of the ith basis element
        w[i] = sum([coeffs[j]*x^j for j in range(len(coeffs))])*xprime/(2*y)  # Replace by local x-coord at infinity
    wint = [(w[i]).integral() for i in range(2*g)]
    return 2*matrix(self.base_ring(), 2*g,2*g, [(w[j]*wint[i]).residue() for i in range(2*g) for j in range(2*g)]) # Residue at \infty_- and \infty_+ is the same, so we can compute only one and multiply by two

def opposite_affine_point(C, P):
    """
    For P=(x,y) on C returns w(P)=(x,-y)
    """
    if P == C(0,1,0):
        return false
    else:
        return C(P[0],-P[1])

def define_alpha_for_infinities(C, prec):
    """
    Returns differential \alpha w.r.t. \omega=x^g/y dx (when one divisor is a difference of two infinities)
    """
    K = C.base_ring()
    p = K.prime()
    f = C.hyperelliptic_polynomials()[0]
    g = C.genus()
    S = monsky_washnitzer.SpecialHyperellipticQuotientRing(C, K)
    MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
    x, y = C.monsky_washnitzer_gens()
    alpha = p * x^(p*g+p-1) * y^(-p) * (1+sum(binomial(QQ(-1)/QQ(2),j)*( (f(x^p)-f(x)^p)*(f(x)^(-p)))**j for j in range(1,prec)))*x.diff()  - p*x^g *y^(-1)*x.diff()
    alpha = MW(alpha)
    return alpha

def reduce_alpha_for_infinities(C, prec):
    """
    Writes alpha=d(pol)+ sum a_i *x^idx/2y where vec=(a_0,...,a_{2g}), where alpha is defined w.r.t. \omega=x^g/y dx (when one divisor is a difference of two infinities)
    """
    alpha = define_alpha_for_infinities(C, prec)
    pol, vec = alpha.reduce_fast()
    return pol, vec

def psi_omega_for_infinities_mixed_basis(C, prec):
    """
    Computes psi(x^gdx/y) in the mixed basis of g holomorphic forms eta_0,...,eta_{g-1} and g Frob^prec remaining ones, i.e., Frob^prec(eta_g),...,Frob^prec(eta_{2g-1})
    """
    K = C.base_ring()
    p = K.prime()
    f = C.hyperelliptic_polynomials()[0]
    g = C.genus()
    basis = h1dr_basis(C)
    pol, vec = reduce_alpha_for_infinities(C, prec)
    M = matrix(2*g,1)
    for i in range(g):
        M[i,0] = vec[i]
    for i in range(g,2*g):
        M[i,0] = vec[i+1]
    Fr = matrix_of_frobenius_on_h1dr(C, p, prec)
    M1 = (Fr-p*identity_matrix(2*g))^(-1) * M
    F = Fr^prec
    F2 = (F.transpose()).matrix_from_rows(range(g,2*g)).list()
    I = identity_matrix(2*g).matrix_from_rows(range(g)).list()
    M2 = ((matrix(2*g,2*g, I + F2)).transpose())^(-1) * M1
    return M2

def height_infinities(C, P, Q, prec):
    """
    Computes h_p(\infty_- - \infty_+, P-Q)=\int_{Q}^P x^gdx/y - sum_{i=0}^{g-1} c_i \int_{Q}^P x^i dx/y, where c_i's are such that psi(x^gdx/y-sum_{i=0}^{g-1} c_i x^i dx/y) is in W = the complementary subspace of holomorphic forms inside H_1^dR(X)
    """
    g = C.genus()
    M = psi_omega_for_infinities_mixed_basis(C, prec)
    PI = C.coleman_integrals_on_basis(Q,P)
    s = 2*PI[g]
    for i in range(g):
        s = s - M[i][0]*PI[i]
    return s 

def tiny_height_infinities(C, P, Q, prec):
    """
    If P and Q are in the same residue disc, computes h_p(\infty_- - \infty_+, P-Q)=\int_{Q}^P x^gdx/y - sum_{i=0}^{g-1} c_i \int_{Q}^P x^i dx/y - integrals are tiny, where c_i's are such that psi(x^gdx/y-sum_{i=0}^{g-1} c_i x^i dx/y) is in W = the complementary subspace of holomorphic forms inside H_1^dR(X)
    """
    g = C.genus()
    M = psi_omega_for_infinities_mixed_basis(C, prec)
    PI = C.tiny_integrals_on_basis(Q,P)
    s = 2*PI[g]
    for i in range(g):
        s = s - M[i][0]*PI[i]
    return s

def height_infinities_withM(C, M, P, Q, prec):
    """
    Computes h_p(\infty_- - \infty_+, P-Q) using the matrix of psi(omega) as an argument M
    """
    g = C.genus()
    PI = C.coleman_integrals_on_basis(Q,P)
    s = 2*PI[g]
    for i in range(g):
        s = s - M[i][0]*PI[i]
    return s

def height_infinities_residuedisc_to_z(C, P, prec):
    """
    Computes h_p(\infty_- - \infty_+, P(z)-P) where P(z) is a point in the residue disc of P
    """
    g = C.genus()
    M = psi_omega_for_infinities_mixed_basis(C, prec)
    PI = C.tiny_integrals_on_basis_to_z(P)
    s = 2*PI[g]
    for i in range(g):
        s = s - M[i][0]*PI[i]
    return s 

def height_infinities_residuedisc_to_z_withM(C, M, P, prec):
    """
    Computes h_p(\infty_- - \infty_+, P(z)-P) where P(z) is a point in the residue disc of P using the matrix of psi(omega) as an argument M
    """
    g = C.genus()
    PI = C.tiny_integrals_on_basis_to_z(P)
    s = 2*PI[g]
    for i in range(g):
        s = s - M[i][0]*PI[i]
    return s 


def height_infinities_first_point_residuedisc_to_z(C, P, Q, prec):
    """
    Computes h_p(\infty_- - \infty_+, P(z)-Q)=h_p(\infty_- - \infty_+, P(z)-P)+h_p(\infty_- - \infty_+, P-Q) where P(z) is a point in the residue disc of P
    """
    return height_infinities_residuedisc_to_z(C, P, prec) + height_infinities(C, P, Q, prec)

def height_infinities_first_point_residuedisc_to_z_withM(C, M, P, Q, prec):
    """
    Computes h_p(\infty_- - \infty_+, P(z)-Q)=h_p(\infty_- - \infty_+, P(z)-P)+h_p(\infty_- - \infty_+, P-Q) where P(z) is a point in the residue disc of P
    """
    return height_infinities_first_point_residuedisc_to_z(C, M, P, prec) + height_infinities_withM(C, M, P, Q, prec)

def Coleman_integral_first_point_residuedisc_to_z(C, P, Q, prec):
    """
    Computes Coleman integrals on basis from fixed Q to P(z), where P(z) is a point in the residue disc of P, uses known integrals from Q to P and tiny integrals on basis from P to P(z)
    """
    return [C.coleman_integrals_on_basis(Q,P)[i] + C.tiny_integrals_on_basis_to_z(P)[i] for i in range(len(C.coleman_integrals_on_basis(Q,P)))]

def height_infinities_nonmonic_changeto_monic(C, P, Q, prec):
    """
    Computes h_p(\infty_- - \infty_+, P-Q) for a non-monic hyperelliptic curve by changing to a monic model
    """
    if C.is_weierstrass(P):
        return 0
    else:
        K = C.base_ring()
        f = C.hyperelliptic_polynomials()[0]
        g = C.genus()
        x = polygen(K)
        X = HyperellipticCurve(sum(f[i]*f[2*g+2]^(2*g+1-i)*x^i for i in range(2*g+3)))
        P1 = X(P[0]*f[2*g+2], P[1]*sqrt(K(f[2*g+2]))^(2*g+1))
        Q1 = X(Q[0]*f[2*g+2], Q[1]*sqrt(K(f[2*g+2]))^(2*g+1))
        return height_infinities(X, P1, Q1, prec)

def psi_forP_H1dR(C, P, prec):
    """
    Computes psi(y(P)/(x-x(P))dx/y) in H^1_{dR} basis for a non-Weierstrass point P
    """
    g = C.genus()
    basis = h1dr_basis(C)
    f = C.hyperelliptic_polynomials()[0]
    if f.degree() == 2*g+1:
        x,y = C.local_coordinates_at_infinity()
        N = C.cup_product_matrix()
    if f.degree() == 2*g+2:
        x,y = local_coordinates_at_infinity_even(C)
        N = cup_product_matrix_even_basis(C, basis)
    wP = P[1] * x.derivative()/((x-P[0])*y)
    NP = opposite_affine_point(C, P)
    PI = C.coleman_integrals_on_basis(NP,P)
    M = matrix(2*g,1)
    w = [0 for i in range(2*g)]
    for i in range(2*g):
        coeffs = basis[i].coeffs()[0][0]  # x-coefficients of the ith basis element
        w[i] = 0
        for j in range(f.degree()-1):
            M[i,0] = M[i,0] + coeffs[j]*PI[j]
            w[i] = w[i] + coeffs[j]*x^j*x.derivative()/(2*y)
        M[i,0] = M[i,0] + (f.degree()-2*g)*(wP*(w[i].integral())).residue()
    return N^(-1)*M 


def psi_forP_mixed_basis(C, P, prec):
    """
    Computes psi(y(P)/(x-x(P))dx/y) in the mixed basis of g holomorphic forms eta_0,...,eta_{g-1} and g Frob^prec remaining ones, i.e., Frob^prec(eta_g),...,Frob^prec(eta_{2g-1})
    """
    g = C.genus()
    K = C.base_ring()
    p = K.prime()
    M = psi_forP_H1dR(C, P, prec)
    f = C.hyperelliptic_polynomials()[0]
    if f.degree() == 2*g+1:
        X = C.change_ring(QQ)
        Fr = X.matrix_of_frobenius(p,prec)
    if f.degree() == 2*g+2:
        Fr = matrix_of_frobenius_on_h1dr(C, p, prec)
    F = Fr^prec
    F2 = (F.transpose()).matrix_from_rows(range(g,2*g)).list()
    I = identity_matrix(2*g).matrix_from_rows(range(g)).list()
    M2 = ((matrix(2*g,2*g, I + F2)).transpose())^(-1) * M    
    return M2
    
def Coleman_Integral_third_kind_antisymmetric(C, P, R, S, prec):
    """
    Computes Coleman integral from S to R of omega_P=y(P)/(x-x(P))dx/y by changing variables to another curve where the integral we integrate becomes x^gdx/y
    """
    if C.is_weierstrass(P):
         return 0
    else:
        K = C.base_ring()
        f = C.hyperelliptic_polynomials()[0]
        g = C.genus()
        x = polygen(K)
        T = (x^(2*g+2) * f(P[0]+1/x) * 1/P[1]^2)
        T.reduce()
        X = HyperellipticCurve(x.parent()(T))
        R1 = X(1/(R[0]-P[0]), (-1)* R[1]/(P[1]*(R[0]-P[0])^(g+1)))
        S1 = X(1/(S[0]-P[0]), (-1)* S[1]/(P[1]*(S[0]-P[0])^(g+1)))
        return 2*X.coleman_integrals_on_basis(S1,R1)[g]

def height_first_antisymmetric(C, P, R, S, prec):
    """
    Computes h_p(P-w(P), R-S) as the integral of omega_P by change of variables and subtract the holomorphic integrals to get the correct differential from the definition of h_p
    """
    g = C.genus()
    IntegralRS = C.coleman_integrals_on_basis(S,R)
    M = psi_forP_mixed_basis(C, P, prec)
    return Coleman_Integral_third_kind_antisymmetric(C, P, R, S, prec) - sum(M[i][0]*IntegralRS[i] for i in range(g))

def height_four_affine_points(C, P, Q, R, S, prec):
    """
    Computes h_p(P-Q, R-S) using and symmetric-antisymmetric decomposition of the first divisor, reduce to the previous case
    """
    K = C.base_ring()
    b = (R[0]-P[0]) * (R[0]-Q[0])^(-1) * (S[0]-Q[0]) * (S[0]-P[0])^(-1)
    logsympart = K(b).log(p_branch=0)
    hasp = height_first_antisymmetric(C, P, R, S, prec)
    hasq = height_first_antisymmetric(C, Q, R, S, prec)
    return 1/2 * (logsympart + hasp - hasq)    

def SplitDivisorsDegreeZero(div):
    """
    Represents div = divisor of degree zero as a sum of (P-Q)
    """
    n = len(div)
    k = sum(div[i][0] for i in range(n))
    if k != 0:
        return "Works only for divisors of degree zero!"
    if n == 2 and div[0][0] == 1 and div[1][0] == -1:
        return [div]
    elif n == 2 and div[0][0] == -1 and div[1][0] == 1:
        return [[div[1], div[0]]]
    else:
        i = 0
        while div[0][0]*div[i][0] > 0:
            i = i + 1
        L = SplitDivisorsDegreeZero([(sgn(div[0][0]),div[0][1]),(sgn(div[i][0]),div[i][1])])
        di = list(div[i])
        di[0] = div[i][0]-sgn(div[i][0])
        if di[0] == 0:
            div.remove(div[i])
        else:
            div[i] = tuple(di)
        d0 = list(div[0])
        d0[0] = div[0][0]-sgn(div[0][0])
        if d0[0] == 0:
            div.remove(div[0])
        else:
            div[0] = tuple(d0)
        return L + SplitDivisorsDegreeZero(div)
    
def height_degreezerodivisors(C, div1, div2, prec):
    """
    Computes h_p(div1, div2) of two divisors of degree zero using the decomposition of each of them
    """
    L1 = SplitDivisorsDegreeZero(div1)
    L2 = SplitDivisorsDegreeZero(div2)
    return sum(height_four_affine_points(C, L1[i][0][1], L1[i][1][1], L2[j][0][1], L2[j][1][1], prec) for i in range(len(L1)) for j in range(len(L2)))

def pair_out(A, B):
    """
    For two lists of the same number of points returns a divisor A-B (of degree zero)
    """
    if len(A) == len(B):
        return [(1, A[i]) for i in range(len(A))] + [(-1, B[j]) for j in range(len(B))]
    else:
        return false

