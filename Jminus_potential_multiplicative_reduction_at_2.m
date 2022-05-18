K<z> := QuadraticField(5);
OK := Integers(K);
p2 := Factorization(2*OK)[1,1];
F<a,c> := FieldOfFractions(PolynomialRing(K,2));
R<x> := PolynomialRing(F);

b := c^5 - a;
f := x^5 - 5*c^2*x^3 + 5*c^4*x - 2*(a-b);

// a has weight 5
// c has weight 1
// the degree of the J-invariants are thus 4, 8, 12, 16, 20

C := HyperellipticCurve(f);
Ig_invs := IgusaInvariants(C);


// Expand in terms of parameter t = a/c^5.
// Assume 2|a, 2 \nmid c.

P<t> := LaurentSeriesRing(K);
Lt := [Evaluate(Ig_invs[i]^5/Ig_invs[5]^i, [t,1]) : i in [1..5]];


// Let Lt[i] = a_n*t^n + a_(n+1)*t^(n+1) + ... for some negative n. It enough to check that the valuation of an is strictly smaller
// from the valuation of a_(n+1) for some i = 1,2,3,4,5.

for i in [1..5] do
	Lti := Lt[i];
	t_val := Valuation(Lti);
	coef_val0 := Valuation(Coefficient(Lti, t_val), p2);
	coef_val1 := Valuation(Coefficient(Lti, t_val+1), p2);
	if (coef_val1 gt coef_val0) and t_val lt 0 then
		printf "Lt%o does always have negative valuation!!!!\n", i;
		break;
	end if;
end for;
