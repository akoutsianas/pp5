/*

In this file, we perform the elimination step for the equation x^p + y^p = z^5 under the assumption 2 \nmid ab and 5 \mid ab with (a,b,c) being a solution of the generalized Fermat equation. In this case, the spaces of Hilbert cusp forms are zero.

*/


Q<z> := PolynomialRing(Rationals());
K<u> := NumberField(z^2 - z - 1);
RK<x> := PolynomialRing(K);
OK := RingOfIntegers(K);
q1 := 1*OK;
q2 := Factorization(2*OK)[1,1];
q5 := Factorization(5*OK)[1,1];


N0 := q1;
N1 := q5;


assert Dimension(NewSubspace(HilbertCuspForms(K, N0))) eq 0;
assert Dimension(NewSubspace(HilbertCuspForms(K, N1))) eq 0;

