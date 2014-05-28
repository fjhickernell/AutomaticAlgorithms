%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
testfun=@(x) sin(2*pi*x);
tau=.1;
epsilon = 1e-2;
a = 0;
b = pi;

%% Kernel
gamma=1;
kernel=@(x,t) exp(-(gamma^2)*(bsxfun(@minus,x,t')).^2);


%% Data and spline approximation
n=15
xnode=linspace(a,b,n)';
Kmat=kernel(xnode,xnode);
condK=cond(Kmat)
y=testfun(xnode);
c=Kmat\y;
splinef=@(x) kernel(x,xnode)*c;

%% Plot of test function and spline
nplot=500;
xplot=linspace(a,b,nplot)';
plot(xplot,testfun(xplot),'b-',...
    xplot,splinef(xplot),'r--',...
    xnode,y,'k.','markersize',20,...
    'linewidth',3)
ntest=1000;
xtest=linspace(a,b,ntest)';

%% Root mean square error
error2=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))
%error_sup=max(abs(testfun(xtest)-splinef(xtest)))
normHsplinef=sqrt(c'*y)
xtPlus=bsxfun(@plus,xnode,xnode');
xtMinus=bsxfun(@minus,xnode,xnode');
%Ktildemat = (sqrt(pi/2)/(2*gamma))*exp(-((gamma^2)/2)*(xtMinus.^2)).*(-erf((gamma/sqrt(2))*(xtPlus-2)) + erf((gamma/sqrt(2))*xtPlus));
Ktildemat = 1/sqrt(2*pi)/(2*gamma)*exp(-gamma^2/2*xtMinus.^2).*(erf(gamma/sqrt(2)*xtPlus)-erf(gamma/sqrt(2)*(xtPlus-2*pi)));
condKtilde=cond(Ktildemat)
%traceK=trace(Kmat)
%traceKt=trace(Ktildemat)
traceKK=trace(Ktildemat/Kmat)
Herrbd=sqrt(1-traceKK)
%guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
%normTsplinef=sqrt(c'*Ktildemat*c);
ErrBound = Herrbd*normHsplinef/(1-tau)
