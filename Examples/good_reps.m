// Find representatives of independent points on Jacobians of hyperelliptic
// curves over Q suitable for p-adic height 
// computations, via the algorithm by Gajovic-Müller.
// This can be used to compute, for instance, the p-adic regulator and
// to solve for the p-adic height pairing in terms of products of
// abelian logs, as required by quadratic Chabauty.
//
// This requires magma version >= 2.25 for MordellWeilGroupGenus2.
//
// For higher genus curves or if your version of magma is older,
// use a subgroup of finite index.


PR<x> := PolynomialRing(RationalField());

function MakeIsZero(R)
  if (Type(R) eq RngSerLaur) then
    if ((R!1) - (R!1) eq 0) then
      return func< z | z eq 0 or RelativePrecision(z) eq 0 >;
    else
      return func< z | z eq 0 or Valuation(z) ge prec > where prec := AbsolutePrecision(R!1);
    end if;
  elif Type(R) eq FldPad then
    return func< z | IsWeaklyZero(z)>;
  else
    return func< z | z eq R!0 >;
  end if;
end function;

function splits_over_Qp(f, Qp)
  fact := Factorisation(ChangeRing(f, Qp));
  if #fact eq 0 then return false; end if;
  return &+[t[2] : t in fact] eq Degree(f);
end function;

function is_pointwise_Qp_rational(P, Qp, Cp)
  // for now assume all pts in L are affine.
  if not splits_over_Qp(P[1], Qp) then
    return false, _;
  end if;
  ap := ChangeRing(P[1], Qp);
  roots := Roots(ap);
  P_list := [];
  for t in roots do
    xx := t[1];
    yy := Evaluate(Parent(ap) ! (P[2]), xx);
    P_list cat:= [[xx, yy] : i in [1..t[2]]];
  end for;
  return true, P_list;
end function;

function bch(A,p)
// Series representation
    //if Valuation(A) gt 10 then A; end if;    
  if IsWeaklyZero(A) then
    S:= "O(";
    S:=S cat IntegerToString(p);
    S:=S cat "^";
    S:=S cat IntegerToString(Precision(Parent(A)));
    S:=S cat ")";
    return S;
  end if;

  v:=Valuation(A);//	if v lt 0 then
  a:=Integers()!(A*p^(-v));//	end if;
  if a lt 0 then
    a+:=p^AbsolutePrecision(A);
  end if;
  b0:=a mod p;
  i:=1;
  bb0:=a;
  L:=[b0];//	b:=(Parent(x)!b0)*x^v;
  while p^i-1 lt Abs(a) do
    bb1:=(bb0-b0) div p;
    b1:= bb1 mod p;
    b0:=b1;
    bb0:=bb1;
    Append(~L,b0);//		b+:=b0*x^(i+v);
    i+:=1;
  end while;
  S:="";
  for i:=1 to #L do
    if L[i] ne 0 then
    S:= S cat IntegerToString(L[i]);
    if i+v-1 ne 0 then
    S:= S cat "*";
    S:= S cat IntegerToString(p);
    S:=S cat "^";
    S:=S cat IntegerToString(i+v-1);
    end if;
    S:=S cat " + ";
    end if;
  end for;
  S:=S cat "O(";
  S:=S cat IntegerToString(p);
  S:=S cat "^";
  S:=S cat IntegerToString(RelativePrecision(A));
  S:=S cat ")";
  return S; //Substring(S,1,#S-2);	
end function;

function sage_points(D)
  // D is a sequence of p-adic points on a hyperelliptic curve C
  // return the points in D in series representation, 
  // suitable for input into Sage.
  p := Prime(Parent(D[1,1]));
  series_list := "[";
  series_list := []; 
  lp := []; 
  for P in D do
    Append(~series_list, " X(" cat bch(P[1], p) cat " , " cat bch(P[2], p) cat " )") ;
  end for;
  return series_list;
end function;

function local_int(C, P, Q, lambda, mu)
  // compute all intersections between divisors representing P and Q
  // See 2014 Math Comp paper by JS Müller
  J := Jacobian(C);
  if P eq Q then 
    lid := LocalIntersectionData(J ! [P[1], P[2]], - J ! [Q[1], Q[2]]: lambda := lambda, mu := mu);
    return [<l[1], -l[2], l[3]> : l in lid];
  else
    return LocalIntersectionData(J ! [P[1], P[2]], J ! [Q[1], Q[2]]: lambda := lambda, mu := mu);
  end if;
end function;

function print_local_int(C, P, Q, lambda, mu, p)
  ints := local_int(C, P, Q, lambda, mu); 
  assert &and[IsOne(t[3]) : t in ints];
  return [[t[1], t[2]] : t in ints | t[1] ne p and t[2] ne 0];
end function;

function rationalize(D)
  rat_div := [ChangeUniverse(Eltseq(P), Rationals()) : P in D];
  return [[P[1], P[2]] : P in rat_div];
end function;

function generators_g2(J)
  // Compute generators of the Mordell-Weil group of a genus 2 Jacobian
  // Uses Stoll's MordellWeilGroupGenus2
  // Requires version >= 2.25.
  A, psi, proved, finite_index := MordellWeilGroupGenus2(J);
  assert proved; // otherwise, more work...
  torsion_orders := [Order(A.i) : i in  [1..#Generators(A)] | Order(A.i) gt 0]; 
  torsion_bas := [psi(A.i) : i in [1..#torsion_orders]];
  free_gens := [psi(A.i) : i in [#torsion_bas+1..#Generators(A)]]; 
  assert #free_gens eq 2; // Need rank = 2
  return torsion_bas, torsion_orders, free_gens;
end function;


function good_representatives(X, p, bas : N := 20, multiple_bound := 10, sage_print := true, curve_name := "");
// X is a hyperelliptic curve with Jacobian of rank g(x)
// Compute representatives of multiples of the points in bas suitable for
// * computing the p-adic regulator of Jac(X)
// * solving for the p-adic height as a bilinear pairing.
// If sage_print is true, the output is also written to a .sage file.
  Fp := GF(p);
  Xred := ChangeRing(X, Fp);
  J := Universe(bas);
  Jred := BaseChange(J, Fp);
  P1red := ProjectiveSpace(Fp,1);
  Qp := pAdicField(p, N);
  XQp := ChangeRing(X, Qp);
  f := HyperellipticPolynomials(X);
  fQp := HyperellipticPolynomials(XQp);
  fred := HyperellipticPolynomials(Xred);
  g := Genus(X);
  r := #bas;

  D_representatives :=  [[] : i in [1..r]];
  E_representatives :=  [[] : i in [1..r]];
  F_representatives :=  [[] : i in [1..r]];
  G_representatives :=  [[] : i in [1..r]];
  lambdas1 := [[0] : i in [1..r]];
  lambdas2 := [[0] : i in [1..r]];
  factors1 := [[0] : i in [1..r]];
  factors2 := [[0] : i in [1..r]];
  
  // First find all good lambdas.
  good_lambdas := [];
  for lambda in [1..p] do
    if IsSquare(Evaluate(fQp, lambda)) and not IsZero(Evaluate(f, lambda)) then
      Append(~good_lambdas, lambda);
    end if;
  end for;
  error if #good_lambdas lt 2, "Not enough good lambdas for this prime";

  multiples := [[n*b : n in [1..multiple_bound]] : b in bas]; // multiples of elements of bas, to avoid recomputing

  bad_factors := []; // n*bas[i] that are bad 

  function is_good(i, n, bad_factors)
    // check if n*bas[i] is bad, meaning that its associated reduced divisor doesn't have pointwise Qp-rational
    // support or its reduction contains a point reducing to a Weierstrass or an infinite disk in its support.
    if <i,n> in bad_factors then 
      return false, bad_factors;
    end if;
    P := multiples[i,n]; 
    if Degree(P[1]) lt g then  // TODO: Get rid of this condition! 
      Append(~bad_factors, <i,n>);
      return false, bad_factors;
    end if;
    b, DP := is_pointwise_Qp_rational(P, Qp, XQp);
    Pred := Jred!P; // TODO: This will cause errors...
    if not b or Degree(Pred[1]) lt g or Degree(GCD(Pred[1], fred)) gt 0  then 
      Append(~bad_factors, <i,n>);
      return false, bad_factors;
    end if;
    // Now check if there's a good lambda to pair with DP
    if &and[IsZero(Evaluate(Pred[1], lambda)) : lambda in good_lambdas] then
      Append(~bad_factors, <i,n>);
      return false, bad_factors;
    end if;
    return true, bad_factors; 
  end function;


  function are_good(P, Q)
    Pred := Jred!P; // TODO: This will cause errors...
    Qred := Jred!Q; // TODO: This will cause errors...
    if Degree(GCD(Pred[1], Qred[1])) gt 0 then 
      // same or opposite disks
      return false, _, _;
    end if;
    // now check if there are good distinct lambdas for both P and Q
    good_lambdas_P := [lambda : lambda in good_lambdas | not IsZero(Evaluate(Pred[1], lambda))];
    good_lambdas_Q := [lambda : lambda in good_lambdas | not IsZero(Evaluate(Qred[1], lambda))];

    function can_choose_distinct(L, M) 
      // L and M are sequences
      // Pick l in L and m in M such that l is not m, if possible
      if #L eq 0 or #M eq 0 then return false, _, _; end if;
      for l in L do 
        for m in M do
          if l ne m then return true, l, m; end if;
        end for;
      end for;
      return false, _, _;
    end function;

    bool, lambda_P, lambda_Q := can_choose_distinct(good_lambdas_P, good_lambdas_Q);
    if not bool then 
      return false, _, _;
    end if;
    return true, lambda_P, lambda_Q;
  end function;

  for i := 1 to r do
    for j := i to r do
      nij := 0;
      mij := 0;
      done := false;
      while not done do 
        nij +:= 1;
        good_nij, bad_factors := is_good(i, nij, bad_factors);
        if good_nij then
          Pij := multiples[i,nij];
          mij := 0;
          while not done and mij lt multiple_bound do 
            mij +:= 1;
            good_mij, bad_factors := is_good(j, mij, bad_factors);
            if good_mij then
              Pji := multiples[j,mij];
              done, lambdaP, lambdaQ := are_good(Pij, Pji);
            end if;
          end while;
        end if;
        error if nij ge multiple_bound, "Need larger multiple_bound";
      end while;          

      factors1[i,j] := nij;
      factors2[i,j] := mij;
      lambdas1[i,j] := lambdaQ; // Note the swap. Since lambdaP is good for P
      lambdas2[i,j] := lambdaP; // we pair DP - E_lambdaQ with DQ - E_lambdaP
      bP, Dij := is_pointwise_Qp_rational(Pij, Qp, XQp);
      bQ, Fij := is_pointwise_Qp_rational(Pji, Qp, XQp);
      D_representatives[i,j] := Dij;
      F_representatives[i,j] := Fij;

      function lambda_reps(lambda)
        blambda, ylambda := IsSquare(Evaluate(fQp, lambda));
        assert blambda;
        Q1 := XQp ! [lambda, ylambda];
        Q2 := XQp ! [lambda, -ylambda];
        Q1_seq := [Eltseq(Q1)[1]/Eltseq(Q1)[3], Eltseq(Q1)[2]/Eltseq(Q1)[3]^(g+1)];
        Q2_seq := [Eltseq(Q2)[1]/Eltseq(Q2)[3], Eltseq(Q2)[2]/Eltseq(Q2)[3]^(g+1)];
        return [Q1_seq, Q2_seq], [Q1, Q2];
      end function;

      E_representatives[i,j] := lambda_reps(lambdas1[i,j]);
      G_representatives[i,j] := lambda_reps(lambdas2[i,j]);

      /* 
       * TODO: If multiple_bound was too small, use this. 
      repeat 
        nij +:= 1;
        if #multiples[i] lt nij then // compute nij*bas[i]
          multiples[i,nij] := bas[i]+multiples[i,nij-1]; 
        end if;
        Dij := multiples[i,nij];
      */
    end for;
  end for;
  
  // Compute intersections for local heights away from p.
  intersections := [[[]] : i in [1..r]];
  away_heights := [[Qp!0: j in [1..r]] : i in [1..r]];
  for i := 1 to r do
    for j := i to r do
      P := multiples[i, factors1[i,j]];
      Q := multiples[j, factors2[i,j]];
      intersections[i,j] := print_local_int(X, P, Q, lambdas1[i,j], lambdas2[i,j], p);
      if #intersections[i,j] ne 0 then 
        away_heights[i,j] := &+[Qp!(t[2])*Log(Qp!t[1]) : t in intersections[i,j]];
      end if;
    end for;
  end for;

  if sage_print then

    SetLogFile(curve_name cat "_p" cat IntegerToString(p) cat "_init.sage" : Overwrite);
    prec1 := Precision(Qp);
    printf "\nR.<x> = Qp(%o,%o)['x']\n", p, prec1;
    printf "X = HyperellipticCurve(%o) \n", HyperellipticPolynomials(X);
    "\ng = X.genus()";
    printf "prechere = %o\n", N; //  TODO!
    printf "p = %o\n", p;
    "K = Qp(p,prechere)";
    "#Mumford coordinates of a QQ-basis of J(QQ):", bas;"";"";
    //"# The following lists [P,Q] are to be read as  divisors P+Q";
    //"\ndef hypinv(pt):\n    return X(pt[0],-pt[1])\n";

    D_series_list := [];
    F_series_list := [];
    E_series_list := [];
    G_series_list := [];
    for i := 1 to r do 
      D_series_list[i] := [sage_points(Dij) : Dij in D_representatives[i]];
      E_series_list[i] := [sage_points(Eij) : Eij in E_representatives[i]];
      F_series_list[i] := [sage_points(Fij) : Fij in F_representatives[i]];
      G_series_list[i] := [sage_points(Gij) : Gij in G_representatives[i]];
      //for j := i to r do 
      // Note the change j -> j-i+1 because D_representatives[i,k] is undefined for k<i
      for j := 1 to r-i+1 do 
        printf "D%o%o = %o\n", i,j+i-1, D_series_list[i,j];
        printf "E%o%o = %o\n", i,j+i-1, E_series_list[i,j];
        printf "F%o%o = %o\n", i,j+i-1, F_series_list[i,j];
        printf "G%o%o = %o\n", i,j+i-1, G_series_list[i,j];
        printf "n%o%o = %o\n", i,j+i-1, factors1[i,j+i-1];
        printf "m%o%o = %o\n", i,j+i-1, factors2[i,j+i-1];
        printf "#Sum of heights away from %o :\naway%o%o = %o\n", p,i,j+i-1, bch(away_heights[i,j+i-1],p);
       // printf "individual intersections :\n %o\n\n\n", intersections[i,j+i-1];
      end for;
    end for;

    UnsetLogFile();
  end if;


  if g eq 2 then 
    // Check that intersections are correct by comparing two different ways to
    // compute the canonical height pairing
    //
          assert HeightPairing(multiples[1,factors1[1,1]],multiples[1,factors2[1,1]]:UseArakelov:=true,lambda:=lambdas1[1,1],mu:=lambdas2[1,1]) 
                - HeightPairing(multiples[1,factors1[1,1]],multiples[1,factors2[1,1]]:UseArakelov:=false) lt 10^-5;
          assert HeightPairing(multiples[1,factors1[1,2]],multiples[2,factors2[1,2]]:UseArakelov:=true,lambda:=lambdas1[1,2],mu:=lambdas2[1,2]) 
                - HeightPairing(multiples[1,factors1[1,2]],multiples[2,factors2[1,2]]:UseArakelov:=false) lt 10^-5;
          assert HeightPairing(multiples[2,factors1[2,2]],multiples[2,factors2[2,2]]:UseArakelov:=true,lambda:=lambdas1[2,2],mu:=lambdas2[2,2])
                - HeightPairing(multiples[2,factors1[2,2]],multiples[2,factors2[2,2]]:UseArakelov:=false) lt 10^-5;
    //assert HeightPairing(-Q,Q:UseArakelov:=true,lambda:=lambda,mu:=mu) - HeightPairing(-Q,Q) lt 10^-5;
  end if;

  return D_representatives, F_representatives, E_representatives, G_representatives, intersections, factors1, factors2;
end function;

