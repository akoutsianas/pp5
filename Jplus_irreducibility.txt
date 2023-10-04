// We prove that the trace of Frobenius at q2 of C^+ is equal to 0


R<x> := PolynomialRing(Rationals());

// In the new model of C^+ after the change of variables (x,y) |-> (1/u, (2*v + 1)/u^3) we divide by 4,
// so we have to consider all classes modulo 8

for a in [1, 3, 5, 7] do
	for c in [0, 2, 4, 6] do
		b := c^5 - a;
		f := (x + 2*c) * (x^5 - 5*c^2*x^3 + 5*c^4*x - 2*(a-b));
		f := (R!(Evaluate(f, 1/x) * x^6) - 1) / 4;
		C := HyperellipticCurve(f, 1);
		J := Jacobian(C);
		Lf := EulerFactor(J, FiniteField(2^2));
		Lf := Reverse(Lf);
		Lfactor := Factorization(R!Lf);
		traceFrob:=[-Coefficient(g[1],1) : g in Lfactor];
		printf "traces: %o\n", traceFrob;
	end for;
end for;
