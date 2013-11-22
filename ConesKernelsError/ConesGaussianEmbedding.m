%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
a=pi;
%testfun=@(x) exp(a*x)-1;
testfun=@(x) (x.^2/2-1/8).*(x>1/2);
tau=2;
epsilon = 1e-2;

%% Kernel
gamma=5;
kernel=@(x,t) exp(-(gamma^2)*(bsxfun(@minus,x,t')).^2);


%% Data and spline approximation
n=5
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
error2=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))
error_sup=max(abs(testfun(xtest)-splinef(xtest)))
%normHsplinef=sqrt(c'*y)
xtPlus=bsxfun(@plus,xnode,xnode');
xtMinus=bsxfun(@minus,xnode,xnode');
Ktildemat = (sqrt(pi/2)/(2*gamma)) ...
    *exp(-((gamma^2)/2)*(xtMinus.^2)) ...
    .*(-erf((gamma/sqrt(2))*(xtPlus-2)) + erf((gamma/sqrt(2))*xtPlus));
condKtilde=cond(Ktildemat)
traceKK=trace(Ktildemat/Kmat)
Herrbd=sqrt(1-traceKK)
%guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
normTsplinef=sqrt(c'*Ktildemat*c);
%htilde=sqrt(2*a^2-trace(Htildemat/Kmat))
ErrBound = tau*Herrbd*normTsplinef/(1-tau*Herrbd)
