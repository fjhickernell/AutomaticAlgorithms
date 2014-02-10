%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
theta=pi;
testfun=@(x) exp(theta*x)-1;
%testfun=@(x) (x.^2/2-1/8).*(x>1/2);
tau=pi;
epsilon = 1e-2;

%% Kernel
a=.5;
b=.9;
kernel=@(x,t) 1-a+a*6*(bp(bsxfun(@plus,acos(x),acos(t'))/(2*pi))+bp(bsxfun(@minus,acos(x),acos(t'))/(2*pi)));

%% Data and spline approximation
n=100
xnode=linspace(1/n,1,n)';
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
error2=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))
error_sup=max(abs(testfun(xtest)-splinef(xtest)))
%normHsplinef=sqrt(c'*y)
xtPlus=bsxfun(@plus,xnode,xnode');
xtTimes=bsxfun(@times,xnode,xnode');
xtMin=bsxfun(@min,xnode,xnode');
xtMax=bsxfun(@max,xnode,xnode');
Ktildemat = gamma^2*(xtMin.^3/3+xtMin.*(xtMax.^2-xtMin.^2)/2+xtTimes.*(1-xtMax));
condKtilde=cond(Ktildemat)
traceKK=trace(Ktildemat/Kmat)
Herrbd=sqrt(gamma/2-traceKK)
%guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
normTsplinef=sqrt(c'*Ktildemat*c);
ErrBound = tau*Herrbd*normTsplinef/(1-tau*Herrbd)
