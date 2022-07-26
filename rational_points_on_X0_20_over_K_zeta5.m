// The modular curve X_0(20).
X := SmallModularCurve(20);

// Base change to L = K(\sqrt{5})).
L := CyclotomicField(5);
EL := ChangeRing(X,L);

// Rank is 0 over L.
assert Rank(EL) eq 0;

// Torsion subgroup over L.
TL, TLtoEL := TorsionSubgroup(EL);

// The map to j-line.
j20 := jInvariant(EL,20);

// Loop through all torsion points over K and check they are cuspidal.
for n := 1 to 6 do
  P := n * TLtoEL(TL.1);
  print P, Evaluate(j20,P);
end for;

