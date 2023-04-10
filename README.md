This code computes the local $p$-adic height $h_p$ on a hyperelliptic curve $C\colon y^2=f(x)$ over $\mathbb{Q}_p$ with good reduction using the construction of Coleman and Gross. More precisely, it computes $h_p(P-Q,R-S)$ for given points $P,Q,R,S$ in $C(\mathbb{Q}_p)$ satisfying the following conditions:

(1) Either all four points are affine, or $\deg(f)$ is even with square leading coefficients and $\{P,Q\} =\{\infty_1,\infty_2\}$ or $\{R,S\} =\{\infty_1,\infty_2\}$, where $\infty_{1/2}\in C(\Q_p)$ are the points at infinity.

(2) The residue discs containing $P$ and $Q$ are distinct from the residue discs containing $R$, $\iota(R)$, $S$, $\iota(S)$, where $\iota$ is the hyperelliptic involution.

(3) Affine points do not lie residue discs at infinity, unless $\{P,Q\}$ or $\{R,S\}$ are in the same residue disc.

The code requires Jennifer Balakrishnan's code for even degree Coleman integrals. To obtain it, download the code from [https://github.com/jbalakrishnan/AWS](https://github.com/jbalakrishnan/AWS) and follow the instructions given there.

Our main commands are the following:

- height_infinities(C, P, Q, prec):    Computes $h_p(\infty_- - \infty_+, P-Q)$

- tiny_height_infinities(C, P, Q, prec):    If $P$ and $Q$ are in the same residue disc, computes $h_p(\infty_- - \infty_+, P-Q)$ up to precision prec using tiny integrals

- height_infinities_withM(C, M, P, Q, prec):     Computes $h_p(\infty_- - \infty_+, P-Q)$ up to precision prec using the matrix of $\psi(\omega)$ (needed for the definition of the local p-adic height) as an argument M, useful if we do a lot of such computations on the same curve

- height_infinities_residuedisc_to_z(C, P, prec):     Computes $h_p(\infty_- - \infty_+, P(z)-P)$ up to precision prec where $P(z)$ is a point in the residue disc of P given parametrically - useful for quadratic Chabauty

- height_infinities_nonmonic_changeto_monic(C, P, Q, prec):    Computes $h_p(\infty_- - \infty_+, P-Q)$ up to precision prec for a non-monic hyperelliptic curve by changing to a monic model if this is possible (i.e., if the leading coefficient is a square), because some algorithms (like reduction in the Monsky-Washnitzer cohomology) are only implemented for monic models of hyperelliptic curves

- height_first_antisymmetric(C, P, R, S, prec):    Computes $h_p(P-\iota(P), R-S)$ up to precision prec as the integral of $\omega_P$ by change of variables in the integral

- height_four_affine_points(C, P, Q, R, S, prec):    Computes $h_p(P-Q, R-S)$ up to precision prec

- height_degreezerodivisors(C, div1, div2, prec):    Computes $h_p(D_1, D_2)$ of two divisors of degree zero using the bi-additivity of the local height via the decomposition into $h_p(P-Q,R-S)$
