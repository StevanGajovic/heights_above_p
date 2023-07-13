This SageMath code computes the local cyclotomic p-adic Coleman--Gross height
pairing $h_p$ on a hyperelliptic curve $C:y^2=f(x)$ over $\mathbb{Q}_p$ with good 
reduction. 
More precisely, it computes $h_p(P-Q,R-S)$ for
$P,Q,R,S\in C(\mathbb{Q}_p)$ satisfying the following conditions:
- $f$ is monic (one can use move_to_monic1 or move_to_monic2 if this is not
  satisfied);
- the residue discs $D(P)$ and $D(Q)$ of $P$ and $Q$ are distinct from
  $D(R)$, $D(\iota(R))$, $D(S)$ and $D(\iota(S))$, where $\iota$ is the
  hyperelliptic involution.

**Dependencies**
Jennifer Balakrishnan's code for even degree Coleman integrals is required.
Download it from [https://github.com/jbalakrishnan/AWS](https://github.com/jbalakrishnan/AWS) and follow the 
instructions given there.

*Main functions*
- height_infinities(P, Q, prec) computes $h_p(\infty_- - \infty_+, P-Q)$
  for even degree models to precision prec. 
- height_four_affine_points(P, Q, R, S, prec) computes $h_p(P-Q, R-S)$ for
  affine $P,Q,R,S$ to precision prec.
- height_divisors(D1, D2) computes $h_p(D_1, D_2)$ for two degree 0
  divisors $D_1,D_2$ on $C$ with disjoint an pointwise $\Q_p$-rational 
- height_infinities_residuedisc_to_z(P) computes $h_p(\infty_- - \infty_+,
  P(z)-P)$, where $P(z)$ is a parametric point in the residue disc of $P$ 
  (used for quadratic Chabauty).
  support.

*Examples*
- $X_0^+(107)$
To compute the rational points on X0+(107), load qc_X0107plus_p7.m into
magma. This requires the QCMod package, available from
https://github.com/steffenmueller/QCMod.
To verify that the coefficients of the global height as a bilinear pairing 
are as claimed in qc_X0107plus_p7.m, run solve_for_height_X0107plus.sage.
- $X_0^+(67)$
To verify the p-adic BSD conjecture for X0+(67) and p = 11, 29, 31, 71, 89,
run bsd_X067plus.sage. 
- Tests 
The file examples.sage contains many tests and sanity checks, including 
comparisons with Balakrishnan's implementation of the algorithm of 
Balakrishnan--Besser for odd degree.

