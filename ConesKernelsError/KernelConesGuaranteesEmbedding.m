%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
testfun=@(x) sin(2*pi*x);
tau=pi;
epsilon = 1e-2;

%% Kernel
a=1;
kernel=@(x,t) a*bsxfun(@min,x,t');

%% Data and spline approximation
n=15;
xnode=linspace(0.01,.99,n)';
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
normHsplinef=sqrt(c'*y)
xtPlus=bsxfun(@plus,xnode,xnode');
xtTimes=bsxfun(@times,xnode,xnode');
xtMin=bsxfun(@min,xnode,xnode');
xtMax=bsxfun(@max,xnode,xnode');
Ktildemat = a^2*(xtMin.^3/3+xtMin.*(xtMax.^2-xtMin.^2)/2+xtTimes.*(1-xtMax));
condKtilde=cond(Ktildemat)
Herrbd=sqrt(a/2-trace(Ktildemat/Kmat))
guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
normTsplinef=sqrt(c'*Ktildemat*c);
%htilde=sqrt(2*a^2-trace(Htildemat/Kmat))
ErrBound = tau*Herrbd*normTsplinef/(1-tau*Herrbd)
