This is ode to compute the local compotent of the p-adic (Coleman-Gross) height above p on a hyperelliptic curve $C:y^2=f(x)$ where p is a prime of good  reduction. For given points $P,Q,R,S$ in $C(\mathbb{Q}_p)$ satisfying the conditions below, computes $h_p(P-Q,R-S)$:
Conditions on points $P,Q,R,S\in C(\mathbb{Q}_p)$:

(1) Either all four points are affine, or one pair of $\{P,Q\}$ and $\{R,S\}$ is the difference $\infty_- -\infty_+$ of points at infinity.

(2) Residue discs of $P$ and $Q$ are distinct of residue discs of $R$, $\iota(R)$, $S$, $\iota(S)$, where $\iota$ is the hyperelliptic involution.

(3) Affine points do not belong to residue discs at infinity, unless $\{P,Q\}$ or $\{R,S\}$ are in the same residue disc.

*** When two points at infinity appear, we consider only the case when $f$ has an even degree, otherwise $f$ can have an arbitrary degree greater than two.

Then (some commands):

height_infinities(C, P, Q, prec):    Computes $h_p(\infty_- - \infty_+, P-Q)$ up to precision prec

tiny_height_infinities(C, P, Q, prec):    If $P$ and $Q$ are in the same residue disc, computes $h_p(\infty_- - \infty_+, P-Q)$ up to precision prec using tiny integrals

height_infinities_withM(C, M, P, Q, prec):     Computes $h_p(\infty_- - \infty_+, P-Q)$ up to precision prec using the matrix of $\psi(\omega)$ (needed for the definition of the local p-adic height) as an argument M, useful if we do a lot of such computations on the same curve

height_infinities_residuedisc_to_z(C, P, prec):     Computes $h_p(\infty_- - \infty_+, P(z)-P)$ up to precision prec where $P(z)$ is a point in the residue disc of P given parametrically - useful for quadratic Chabauty

height_infinities_nonmonic_changeto_monic(C, P, Q, prec):    Computes $h_p(\infty_- - \infty_+, P-Q)$ up to precision prec for a non-monic hyperelliptic curve by changing to a monic model if this is possible (i.e., if the leading coefficient is a square), because some algorithms (like reduction in the Monsky-Washnitzer cohomology) are only implemented for monic models of hyperelliptic curves

height_first_antisymmetric(C, P, R, S, prec):    Computes $h_p(P-\iota(P), R-S)$ up to precision prec as the integral of $\omega_P$ by change of variables in the integral

height_four_affine_points(C, P, Q, R, S, prec):    Computes $h_p(P-Q, R-S)$ up to precision prec

height_degreezerodivisors(C, div1, div2, prec):    Computes $h_p(D_1, D_2)$ of two divisors of degree zero using the bi-additivity of the local height via the decomposition into $h_p(P-Q,R-S)$
