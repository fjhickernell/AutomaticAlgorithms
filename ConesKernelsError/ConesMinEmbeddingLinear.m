%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
testfun = @(x) sin(2*pi*x);
theta=pi;
tau=.1;
epsilon = 1e-2;


%% Kernel
gamma=1;
kernel=@(x,t) gamma*bsxfun(@min,x,t');
%% Data and spline approximation
n=20
xnode=linspace(-1,1,n)';
Kmat=kernel(xnode,xnode);
condK=cond(Kmat)
y=testfun(xnode);
c=Kmat\y;
splinef=@(x) kernel(x,xnode)*c;

%% Plot of test function and spline
ntest=1000;
xtest=linspace(-1,1,ntest)';
ytest = interp1q(xnode,y,xtest);
plot(xtest,ytest,'b-',...
    xtest,splinef(xtest),'r--',...
    xnode,y,'k.','markersize',20,...
    'linewidth',3)

%% Root mean square error
error2=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))
error_sup=max(abs(testfun(xtest)-splinef(xtest)))
normHsplinef=sqrt(c'*y)
xtPlus=bsxfun(@plus,xnode,xnode');
xtTimes=bsxfun(@times,xnode,xnode');
xtMin=bsxfun(@min,xnode,xnode');
xtMax=bsxfun(@max,xnode,xnode');
Ktildemat = gamma^2*(xtMin.^3/3+xtMin.*(xtMax.^2-xtMin.^2)/2+xtTimes.*(1-xtMax));
condKtilde=cond(Ktildemat)
traceKK=trace(Ktildemat/Kmat)
Herrbd=sqrt(gamma/2-traceKK)

%% Algorithm 1 Stage 1
%normTsplinef=sqrt(c'*Ktildemat*c);
ErrBound = Herrbd*normHsplinef/(1-tau)
