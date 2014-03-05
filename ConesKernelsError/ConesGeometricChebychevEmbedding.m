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
a=.5;
b=.9;
kernel=@(x,t) 1-a+a*(1-b)*((bsxfun(@times,x,t')-bsxfun(@times,sqrt(1-x.^2),sqrt(1-t'.^2))-b)./(1+b^2-2*b*bsxfun(@times,x,t')-bsxfun(@times,sqrt(1-x.^2),sqrt(1-t'.^2))) ...
    +(bsxfun(@times,x,t')+bsxfun(@times,sqrt(1-x.^2),sqrt(1-t'.^2))-b)./(1+b^2-2*b*bsxfun(@times,x,t')+bsxfun(@times,sqrt(1-x.^2),sqrt(1-t'.^2))));

%% Data and spline approximation
n=100
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
%normHsplinef=sqrt(c'*y)
xtTimes = bsxfun(@times,xnode,xnode');
xtsTimes = bsxfun(@times,sqrt(1-xnode.^2),sqrt(1-xnode'.^2));
xtMxts = xtTimes-xtsTimes;
xtPxts = xtTimes+xtsTimes;
Ktildemat = (1-a)^2+a^2*(1-b)^2*((xtMxts-b^2)./(1+b^4-2*b^2*xtMxts)+(xtPxts-b^2)/(1+b^4-2*b^2*xtPxts));
condKtilde=cond(Ktildemat)
traceKK=trace(Ktildemat/Kmat)
Herrbd=sqrt(1-traceKK)
%guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
normTsplinef=sqrt(c'*Ktildemat*c);
ErrBound = tau*Herrbd*normTsplinef/(1-tau*Herrbd)
