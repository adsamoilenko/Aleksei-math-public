p=7;
t(p)=exp(2*Pi*I/p);
s(p)=exp(2*Pi*I/(p-1));

B1(x)=t(p)^(lift(x));
B2(x)=sum(z=1,p-1,B1(Mod(x/z,p))*B1(z));
DUPL(x)=B2(x)^2;
MULT(x,y)=B2(x)*B2(y);

Fourier(f,y,p)=sum(i=1,p-1,f(i)*s(p)^(znlog(i,znprimroot(p))*y));
twoFourier(f,x,y,p)=sum(i=1,p-1,sum(j=1,p-1,f(i,j)*s(p)^(znlog(i,znprimroot(p))*x)*s(p)^(znlog(j,znprimroot(p))*y)));
invFourier(F,x,p)=(1/(p-1))*sum(i=1,p-1,F(i)*s(p)^(-znlog(x,znprimroot(p))*i));
invtwoFourier(F,u,v,p)=(1/(p-1)^2)*sum(i=1,p-1,sum(j=1,p-1,F(i,j)*s(p)^(-znlog(u,znprimroot(p))*i)*s(p)^(-znlog(v,znprimroot(p))*j)));

hatB2(y)=Fourier(B2,y,p);
hatDUPL(y)=Fourier(DUPL,y,p);
hatDKER(y)=hatDUPL(y)/hatB2(y);
DKER(x)=invFourier(hatDKER,x,p);
for(i=1,p-1,print(i, " | ", bestappr(DKER(i)), " | ", p+kronecker(1-4*i,p)));

hatMULT(x,y)=twoFourier(MULT,x,y,p);
diagB2(x,y)=if(x==y,B2(x),0);
hatdiagB2(x,y)=twoFourier(diagB2,x,y,p);
hatMKER(x,y)=hatMULT(x,y)/hatdiagB2(x,y);
MKER(x,y)=invtwoFourier(hatMKER,x,y,p);
forvec(I=[p-1,p-1],X=I[1]+1;Y=I[2]+1;print(X, Y ," | ", bestappr(MKER(X,Y)), " | ", -1+kronecker((X)^2+(Y)^2+1-2*(X*Y+X+Y),p);));
