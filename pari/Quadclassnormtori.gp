Quadclassnormtori(l,m,n)=
{
local(kn,Lnkn,Ln,UULn,M,epsn,taun);
kn=bnfinit(polcyclo(l^n,y));
Lnkn=rnfinit(kn,x^2+m,1);
Ln=bnfinit(Lnkn[11][1]);
my(hLn = Ln.no);
my(hkn = kn.no);

UULn=bnfunits(Ln);
M=Mat([]);
for(i=1,Ln.r1+Ln.r2,M=concat(M,bnfisunit(kn,rnfeltnorm(Lnkn,nffactorback(Ln,UULn[1][i])))));
M*Mod(1,2);
epsn=kn.r2-matrank(M*Mod(1,2));

taun=#idealfactor(kn,quaddisc(-m));

my(hn  = hLn/(hkn*2^(taun-epsn-1)));

return([["eps(n)","tau(n)","h(Ln)","h(kn)","h(n)"]; [epsn,taun,hLn,hkn,hn]]);
}