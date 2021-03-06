%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
theta=pi;
%testfun=@(x) exp(theta*x)-1;
testfun=@(x) exp(-5*x).*sin(10*x);
tau=2;
epsilon = 1e-2;

%% Kernel
n=20
a=2;
b=10^(-3/n);
kernel=@(x,t) a*sqrt((1-b)/(1+b))*exp(-a^2*b*(bsxfun(@minus,x,t')).^2/(1-b^2)+1/2*(1-a^2*(1-b)/(1+b))*bsxfun(@plus,x.^2,t'.^2));

%% Data and spline approximation
xnode=linspace(-1,1,n)';
xnode=sin(xnode*pi/2);
Kmat=kernel(xnode,xnode);
condK=cond(Kmat)
y=testfun(xnode);
c=Kmat\y;
splinef=@(x) kernel(x,xnode)*c;

%% Plot of test function and spline
nplot=500;
xplot=linspace(-1,1,nplot)';
plot(xplot,testfun(xplot),'b-',...
    xplot,splinef(xplot),'r--',...
    xnode,y,'k.','markersize',20,...
    'linewidth',3)
ntest=1000;
xtest=linspace(-1,1,ntest)';

%% Root mean square error
error2=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))
error_sup=max(abs(testfun(xtest)-splinef(xtest)))
%normHsplinef=sqrt(c'*y)
Ktildemat = a*(1-b)^2/sqrt(1-b^4)*exp(-a^2*b^2*(bsxfun(@minus,xnode,xnode')).^2/(1-b^4)+1/2*(1-a^2*(1-b^2)/(1+b^2))*bsxfun(@plus,xnode.^2,xnode'.^2));
condKtilde=cond(Ktildemat)
traceKK=trace(Ktildemat/Kmat)
Herrbd=sqrt(1-traceKK)
%guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
normTsplinef=sqrt(c'*Ktildemat*c);
ErrBound = tau*Herrbd*normTsplinef/(1-tau*Herrbd)
