%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
theta=pi;
testfun=@(x) exp(theta*x)-1;
%testfun=@(x) sin(2*pi*x);
%testfun=@(x) (x.^2/2-1/8).*(x>1/2);
tau=.1;
epsilon = 1e-2;

%% Kernel
a=.5;
%b=.9;
bp2 = @(x) x.^2-x+1/6;
bp4 = @(x) x.^4-2*x.^3+x.^2-1/30;
kernel=@(x,t) 1-a+a*6*(bp2(abs(bsxfun(@plus,acos(x),acos(t')))/(2*pi))+bp2(abs(bsxfun(@minus,acos(x),acos(t')))/(2*pi)));

%% Data and spline approximation
n=80
xnode=linspace(-1,1,n)';
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
normHsplinef=sqrt(c'*y)
Ktildemat = (1-a)^2-12*a^2*(bp4(abs(bsxfun(@plus,acos(xnode),acos(xnode')))/(2*pi))+bp4(abs(bsxfun(@minus,acos(xnode),acos(xnode')))/(2*pi)));
condKtilde=cond(Ktildemat)
traceKK=trace(Ktildemat/Kmat)
Herrbd=sqrt(1-traceKK)

%% Algorithm 1 Stage 1
%normTsplinef=sqrt(c'*Ktildemat*c);
ErrBound = Herrbd*normHsplinef/(1-tau)
