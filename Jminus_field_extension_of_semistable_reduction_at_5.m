// Confirm the ramification degree of the extension over Jminus get semistable reduction at 5

Q5 := pAdicField(5,100);
RQ<u> := PolynomialRing(Q5);
K<z> := TotallyRamifiedExtension(Q5,u^2-5);
OK := Integers(K);
RK<s> := PolynomialRing(K);
M := TotallyRamifiedExtension(K,s^4-z);
OM := Integers(M);
RM<t> := PolynomialRing(M);


N:=25;

for a0:=1 to N do
	for c0:=1 to N do
		b0 := c0^5 - a0;
		d0 := -4*a0 - a0^5 - 10*a0^4*c0 - 35*a0^3*c0^2 - 50*a0^2*c0^3 - 25*a0*c0^4;

		if a0*b0 mod 5 ne 0 then
			if Valuation(d0,5) eq 1 then
    				print "-> L case";
				L := RamifiedRepresentation(LocalField(M,t^5-5*c0^2*t^3+5*c0^4*t-2*(a0-b0)));
				OL := Integers(L);
				print RamificationDegree(L) * RamificationDegree(M);
			elif Valuation(d0,5) ge 2 then
				print "-> M case";
				print RamificationDegree(M);
			elif Valuation(d0,5) eq 0 then
				assert false;
			end if;
		end if;
	end for;
end for;


