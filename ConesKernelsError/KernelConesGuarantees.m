%Cones and Kernel Approximation for Guarantees
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
testfun=@(x) sin(8*x);

a=10;
kernel=@(x,t) exp(-(a.^2)*(bsxfun(@minus,x,t')).^2);

n=20;
xnode=linspace(0,1,n)';
Kmat=kernel(xnode,xnode);
%xnode=0.5;
y=testfun(xnode);
c=Kmat\y;
splinef=@(x) kernel(x,xnode)*c;

nplot=500;
xplot=linspace(0,1,nplot)';
plot(xplot,testfun(xplot),'b-',...
    xplot,splinef(xplot),'r--',...
    xnode,y,'k.','markersize',20,...
    'linewidth',3)
ntest=1000;
xtest=linspace(0,1,ntest)';
error=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))

normHsplinef=c'*y
xpxt=bsxfun(@plus,xnode,xnode');
invKmat=inv(Kmat);
Ktildemat=(sqrt(pi/2)/(2*a)) ...
    *exp(-((a^2)/2)*(bsxfun(@minus,xnode,xnode').^2)) ...
    .*(-erf((a/sqrt(2))*(xpxt-2)) + erf((a/sqrt(2))*xpxt));
Herrbd=1-sum(sum(Ktildemat'.*invKmat,2))
guesserrest=Herrbd*normHsplinef
