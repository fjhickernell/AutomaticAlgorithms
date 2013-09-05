%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
testfun=@(x) sin(8*x);

%% Kernel
a=20;
kernel=@(x,t) exp(-a*abs(bsxfun(@minus,x,t'))).*(1+abs(bsxfun(@minus,x,t')));

%% Data and spline approximation
n=25;
xnode=linspace(0,1,n)';
Kmat=kernel(xnode,xnode);
condK=cond(Kmat)
y=testfun(xnode);
c=Kmat\y;
splinef=@(x) kernel(x,xnode)*c;

%% Plot of test function and spline
nplot=500;
xplot=linspace(0,1,nplot)';
plot(xplot,testfun(xplot),'b-',...
    xplot,splinef(xplot),'r--',...
    xnode,y,'k.','markersize',20,...
    'linewidth',3)
ntest=1000;
xtest=linspace(0,1,ntest)';

%% Root mean square error
error=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))

normHsplinef=c'*y
xpxt=bsxfun(@plus,xnode,xnode');
Ktildemat=(sqrt(pi/2)/(2*a)) ...
    *exp(-((a^2)/2)*(bsxfun(@minus,xnode,xnode').^2)) ...
    .*(-erf((a/sqrt(2))*(xpxt-2)) + erf((a/sqrt(2))*xpxt));
condKtilde=cond(Ktildemat)
Herrbd=1-trace(Ktildemat/Kmat)
guesserrest=Herrbd*normHsplinef
