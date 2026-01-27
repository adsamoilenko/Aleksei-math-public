p=37;
k=1;
F=bnfinit(polcyclo(p));
zz=Mod(x,polcyclo(p));

S=idealprimedec(F,p);
US=bnfunits(F,S);

uu_chi(d)=prod(a=1,p-1,(1-zz^a)^(lift(Mod(1/(p-1),p^k))*lift(teichmuller(a+O(p^k))^(-d))));

for(d=1,p-1,print(d," -valuations of exponents ", vv=vector(length(bnfisunit(F,uu_chi(d),US)), i, valuation(bnfisunit(F,uu_chi(d),US)[i], p))," -min- ",vecmin(vv)));
