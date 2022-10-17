
K := QuadraticField(5);
OK := Integers(K);
I2 := Factorization(2*OK)[1,1];
R<x> := PolynomialRing(K);
d := 1;
M := 2^6;


// Case d = 1, a mod 32 = 0

for a:=1 to M do
	for b:=1 to M do

		if a mod 32 eq 0 and b mod 4 eq 3 then
			E := QuadraticTwist(EllipticCurve(x*(x+a)*(x-b)),d);
			LI := LocalInformation(E,I2);
			assert LI[3] eq 1;
		end if;
	end for;
end for;

// Case d = 1, b mod 32 = 0

for a:=1 to M do
	for b:=1 to M do
		if b mod 32 eq 0 and a mod 4 eq 1 then
			E := QuadraticTwist(EllipticCurve(x*(x+a)*(x-b)),d);
			LI := LocalInformation(E,I2);
			assert LI[3] eq 1;
		end if;
	end for;
end for;

