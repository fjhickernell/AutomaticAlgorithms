clear all, close all
sigma=0.001;
t=rand;
f=@(x) exp(-(x-t).^2/(2*sigma^2))...
    /(sqrt(2*pi)*sigma);
chebf=chebfun(f,[0 1]);

%Exact integral
tic
exactint=normcdf((1-t)/sigma)-normcdf(-t/sigma)
toc
tol=1e-14;

%Approximate integrals
tic
appxint=quad(f,0,1,tol) %quad
toc
err=exactint-appxint

tic
appxintgk=quadgk(f,0,1,'AbsTol',tol) %quadgk
toc
errgk=exactint-appxintgk

tic
appxintcheb=sum(chebf,0,1) %chebfun
toc
errcheb=exactint-appxintcheb

%Function approximation
nxx=1e5;
xx=rand(nxx,1);
fxx=f(xx);

ndat=1000;
xdat=(0:ndat-1)/(ndat-1);
ydat=f(xdat);
tic
yspline=spline(xdat,ydat,xx);
toc
err=abs(fxx-yspline);
maxerrspline=max(err)
rmserrspline=sqrt(mean(err.^2))

tic
ycheb=chebf(xx);
toc
err=abs(fxx-ycheb);
maxerrcheb=max(err)
rmserrcheb=sqrt(mean(err.^2))



