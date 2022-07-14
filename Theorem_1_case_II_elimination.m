/*

In this file, we perform the elimination step for the equation x^p + y^p = z^5 under the assumption 2 | ab and 5 \nmid ab with (a,b,c) being a solution of the generalized Fermat equation. 

*/


Q<z> := PolynomialRing(Rationals());
K<u> := NumberField(z^2 - z - 1);
r5 := 2*u - 1;
RK<x> := PolynomialRing(K);
OK := RingOfIntegers(K);
q2 := SetToSequence(Support(2*OK))[1];
q5 := SetToSequence(Support(5*OK))[1];
Ds := [1, K!(-1), K!2, K!(-2), -2 + 2*r5, -2 - 2*r5, -4 + 4*r5, -4 - 4*r5];

// Q5<z5> := ext<K | x^2 + u + 2>;
// OQ5 := RingOfIntegers(Q5);


N0 := q2 * q5^2;
N1 := q2 * q5^3;


decomp0 := NewformDecomposition(NewSubspace(HilbertCuspForms(K, N0)));
decomp1 := NewformDecomposition(NewSubspace(HilbertCuspForms(K, N1)));


function get_traceFrob(C)
	// The traces of Frobenius with respect to the primes above the rational prime q
	
	Zq := ZetaFunction(C);
	Pq := Numerator(Zq);
	P := Numerator(Evaluate(Pq,1/z)*z^4);
	traceFrob := [];
	for fac in Factorization(ChangeRing(P,K)) do
		Append(~traceFrob, -Coefficient(fac[1],1));
	end for;

	return traceFrob;
end function;


function eliminate_newform(fnew : Qs := [3, 7, 11])
	// We apply elimination step for a given Hilbert newform fnew

	small_primes := [];
	for d in Ds do
		Td := [];
		for q in Qs do
			Tq := q;
			qK := SetToSequence(Support(q*OK))[1];
			aqKfnew := HeckeEigenvalue(Eigenform(fnew), qK);
			bol, KfnewtoK := IsSubfield(HeckeEigenvalueField(fnew), K);
			assert bol;
			Fq, OKtoFq := ResidueClassField(qK);
			_<t> := PolynomialRing(Fq);
			dbar := OKtoFq(d);
			for a, c in [1..q] do
				b := c^5 - a;
				if (a*b) mod q ne 0 then
					C := HyperellipticCurve(t^5 - 5*c^2*t^3*dbar^2 + 5*c^4*t*dbar^4 - 2*(a-b)*dbar^5);
					traceFrob := get_traceFrob(C);
					Tq *:= Norm(&*[KfnewtoK(aqKfnew) - aqKJ: aqKJ in traceFrob]);
				else
					Tq *:= Norm((Norm(qK) + 1)^2 - aqKfnew^2);
				end if;
			end for;
			Append(~Td, Integers()!Tq);
		end for;

		for q in PrimeFactors(Gcd(Td)) do
			if q notin small_primes then
				Append(~small_primes, q);
			end if;
		end for;
	end for;

	return small_primes;
end function;



// We apply the elimination step for the Level q2*q5^2

for newf in decomp0 do
	newf_small_primes := eliminate_newform(newf);
	printf "newf_small_primes: %o\n", newf_small_primes;
end for;


// We apply the elimination step for the Level q2*q5^3

for newf in decomp1 do
        newf_small_primes := eliminate_newform(newf);
        printf "newf_small_primes: %o\n", newf_small_primes;
end for;


