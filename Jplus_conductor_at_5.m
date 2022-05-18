K<r5> := QuadraticField(5);
F<a,b,c> := FieldOfFractions(PolynomialRing(K,3));
R<x> := PolynomialRing(F);
OK := RingOfIntegers(K);
p5 := ideal<OK | r5>;
Fp := FiniteField(5);
_<t> := PolynomialRing(Fp);

f := (x + 2*c) * (x^5 - 5*c^2*x^3 + 5*c^4*x - 2*(a-b));
fa := Evaluate(f, r5*x - 2*c - a) / r5^6;
fb := Evaluate(f, r5*x - 2*c + a) / r5^6;

// We compute the valuation vectors and check double root criterion for all possible classes modulo 5 of (a, b, c)

// Case 5|a

printf "Case 5|a";
for ai in [0..4] do
	for ci in [1..4] do
		coefs := [Evaluate(Coefficient(fa, k), [5^5*ai, ci^5 - 5^5*ai, ci]) : k in [0..6]];
		val_vecs := [Valuation(coefs[k], p5) : k in [1..7]];
		printf "valuation vectors: %o\n", val_vecs;
		bfa := &+[Fp!(OK!coefs[k] mod p5) * t^(k-1): k in [1..7]];
		printf "Double root is satisfied: %o\n", Max([rt[2] : rt in Roots(bfa)]) eq 2;
	end for;
end for;


// Case 5|b
printf "Case 5|b";
for ai in [0..4] do
        for ci in [1..4] do
		coefs := [Evaluate(Coefficient(fb, k), [5^5*ai, ci^5 - 5^5*ai, ci]) : k in [0..6]];
                val_vecs := [Valuation(coefs[k], p5) : k in [1..7]];
                printf "valuation vectors: %o\n", val_vecs;
                bfa := &+[Fp!(OK!coefs[k] mod p5) * t^(k-1): k in [1..7]];
                printf "Double root is satisfied: %o\n", Max([rt[2] : rt in Roots(bfa)]) eq 2;
        end for;
end for;

