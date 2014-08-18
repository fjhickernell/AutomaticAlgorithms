%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
theta=pi;
testfun=@(x) exp(theta*x)-1;
%testfun=@(x) (x.^2/2-1/8).*(x>1/2);
tau=2;
epsilon = 1e-2;

%% Kernel
a = .5;
b = 2;
bp2 = @(x) x.^2-x+1/6;
bp4 = @(x) x.^4-2*x.^3+x.^2-1/30;
bp8 = @(x) x.^8-4*x.^7+14*x.^6/3-7*x.^4/3+2*x.^2/3-1/30;
kernel=@(x,t) 1-a+a*(-1)^(b+1)*(2*pi)^(2*b)/(2*factorial(2*b)*zeta(2*b))*(bp4(abs(bsxfun(@plus,acos(x),acos(t')))/(2*pi))+bp4(abs(bsxfun(@minus,acos(x),acos(t')))/(2*pi)));

%% Data and spline approximation
n=40
%xnode=linspace(-1,1,n)';
xnode = cos(linspace(pi/(2*n),(2*n-1)/(2*n)*pi,n)');
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
Ktildemat = (1-a)^2-a^2*(2*pi)^(4*b)/(2*factorial(4*b)*zeta(2*b)^2)*(bp8(abs(bsxfun(@plus,acos(xnode),acos(xnode')))/(2*pi))+bp8(abs(bsxfun(@minus,acos(xnode),acos(xnode')))/(2*pi)));
condKtilde=cond(Ktildemat)
traceKK=trace(Ktildemat/Kmat)
Herrbd=sqrt(1-traceKK)
%guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
normTsplinef=sqrt(c'*Ktildemat*c)
ErrBound = tau*Herrbd*normTsplinef/(1-tau*Herrbd)
