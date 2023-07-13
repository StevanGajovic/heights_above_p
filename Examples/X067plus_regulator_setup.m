// Compute good representatives to compute p-adic regulators
// on J0(67)+
//

load "good_reps.m";
C := HyperellipticCurve(x^6 + 4*x^5 + 2*x^4 + 2*x^3 + x^2 - 2*x + 1);


torsion_bas, torsion_orders, free_gens := generators_g2(Jacobian(C));
assert #torsion_bas eq 0;

Dreps := [**];
Ereps := [**];
Freps := [**];
Greps := [**];
ints := [**];
facts1 := [**];
facts2 := [**];


good_even_ordinary_primes := [ 11, 29, 31, 71, 89];
// primes p below 100 such that
// - J0(67)+ has good ordinary reduction at p
// - Sage can compute the p-adic L function, i.e.
//   - p splits in the RM field Q(sqrt(5))
//   - the modular symbol used in the p-adic L-function computation
//     is defined over Qp
// - X0(67)+(Qp) contains no Weierstrass point
//
D_representatives := [];
E_representatives := [];
F_representatives := [];
G_representatives := [];
intersections := [];
factors1 := [];
factors2 := [];

// Let free_gens = [P1, P2]. This computes divisors of degree 2 and
// integers such that
// * [D11-E11] = n11*P1, [F11-G11] = m11*P1, 
// * [D12-E12] = n12*P1, [F12-G12] = m12*P2, 
// * [D22-E22] = n22*P1, [F22-G22] = m22*P2, 
// * all divisors are supported in affine points
// * the E-divisors and G-divisors are canonical.
// It also computes the sum of all local heights away from p:
//  h_v(Dij-Eij, Fij-Gij).
//
// Besides returning these data, they are also written to a .sage file.
// 


for p in good_even_ordinary_primes do
  try
      
    D_representatives, F_representatives, E_representatives, G_representatives, 
      intersections, factors1, factors2 := 
  good_representatives(C, p, free_gens : multiple_bound := 20, sage_print, curve_name := "X067plus");

  catch e;
    e;
  end try;
end for;
