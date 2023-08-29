//
// Compute data to determine the 7-adic height as a bilinear pairing on 
// J0(107)+. In particular, find suitable representatives to which we can 
// apply Algorithm 3.1 in arXiv:2307.15787. Write this data to a .sage file.

load "good_reps.m";
_<x> := PolynomialRing(Rationals());
 

f107 := x^6 + 2*x^5 + 5*x^4 + 2*x^3 - 2*x^2 - 4*x - 3;
X := HyperellipticCurve(f107); 

J := Jacobian(X);

N := 12;
f := HyperellipticPolynomials(X); 
gX := Genus(X);
ptsX := Points(X:Bound:=100);
"Small points: ", ptsX;


// Compute generators for the Mordell-Weil group 
torsion_bas, torsion_orders, bas := generators_g2(J);
printf "The rank is %o.\n", #bas; // rank = 2
// This spares us the trouble of checking saturation in MW sieve computation.
bas[2] := -bas[2]; // These generators work better for this example.


p := 7;
D_representatives, F_representatives, E_representatives, G_representatives, intersections, factors1, factors2 := 
         good_representatives(X, p, bas : N := N, multiple_bound := 20,  sage_print := true, curve_name := "X0107plus");


