// This program checks that the image of units of OK under reduction mod \Fq_2^2 is surjective.

K<z> := QuadraticField(5);
OK := Integers(K);
I2 := Factorization(2*OK)[1,1];

UK,f := UnitGroup(OK);
u := K!f(UK.2);

R2 := quo<OK | 4>;
U2, f2 := UnitGroup(R2);

powers := {};
for n:=-10 to 10 do
	powers := powers join {(u^n)@@f2} join {(-u^n)@@f2};
end for;

assert #powers eq #U2;
