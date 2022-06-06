// The conductor computations of J^- at 5

Q5 := pAdicField(5,100);
RQ<u> := PolynomialRing(Q5);
K<r5> := TotallyRamifiedExtension(Q5, u^2-5);
RK<s> := PolynomialRing(K);
M := TotallyRamifiedExtension(K, s^4-r5);
RM<t> := PolynomialRing(M);


// The conductor when C^- attains good reduction over the 4 degree extension M/K

G := GaloisRepresentations(M,K);
for rho in G do
	if IsFaithful(Character(rho)) then
		assert Valuation(Conductor(rho)) eq 1;
	end if;
end for;


// The conductor when C^- attains good reduction over the 20 degree extension L/K

N:=25;
for a0:=1 to N do
	for c0:=1 to N do
		b0 := c0^5 - a0;
		d0 := -4*a0 - a0^5 - 10*a0^4*c0 - 35*a0^3*c0^2 - 50*a0^2*c0^3 - 25*a0*c0^4;
		if (a0*b0 mod 5 ne 0) and (Valuation(d0,5) eq 1 ) then
			print "*****";
			printf "a0:%o, b0:%o, c0:%o\n", a0 mod 5,b0 mod 5,c0 mod 5;
			print "-> L case";
			Q5 := pAdicField(5,100);
			L := RamifiedRepresentation(LocalField(M, t^5-5*c0^2*t^3+5*c0^4*t-2*(a0-b0)));
			RL<x> := PolynomialRing(L);
			OL := Integers(L);
			OK := Integers(K);
			G := GaloisRepresentations(L,K);
    			for rho in G do
      				if IsFaithful(Character(rho)) then
       					assert Valuation(Conductor(rho)) eq 3;
				end if;
			end for;
    			print "The conductor exponent is = 3.";
		end if;
	end for;
end for;
