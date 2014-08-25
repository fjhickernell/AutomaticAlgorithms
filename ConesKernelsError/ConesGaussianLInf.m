%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
n = 100
a = 0; b = 1;
xnode=linspace(a,b,n)';
h = (b-a)/(2*n); R = (b-a)/2;
testfun=@(x) sin(pi*x);
gamma = 1;
kernel=@(x,t) exp(-gamma*(bsxfun(@minus,x,t')).^2);
tau = 2;
%c0 = R/6; C1 = 2304*gamma; C2 = log(4)+16; C3 = C1*exp(C2);
%C = min(c0/2,1/(exp(1)*C3))/2;
%C = min(c0/2,1/sqrt(exp(1)*C3))/2;
C4 = 64/3; %C5 = 1/2; c1 = sqrt(2*gamma^2*C5*C4^3*(6+C5*C4));
C = 3*sqrt(2)*gamma*C4^2;
h0 = min([1/sqrt(tau*C) R/C4 sqrt(3)/C4/sqrt(2*gamma)]);
%epsilon = 1e-2;

%% Data and spline approximation
Kmat=kernel(xnode,xnode);
condK=cond(Kmat)
y=testfun(xnode);
c=Kmat\y;
splinef=@(x) kernel(x,xnode)*c;

%% Plot of test function and spline
ntest=5000;
xtest=linspace(a,b,ntest)';
ytest = interp1q(xnode,y,xtest);
plot(xtest,ytest,'b-',...
    xtest,splinef(xtest),'r--',...
    xnode,y,'k.','markersize',20,...
    'linewidth',3)

%% Root mean square error
%error2=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))
error_sup=max(abs(testfun(xtest)-splinef(xtest)))
normHsplinef=sqrt(c'*y)


%% Algorithm 1 Stage 1
%[~,fval] = fminbnd(@(x) -FcnToMax(x,xnode,c,gamma,1),a,b);
%eBoundEst = -tau*exp(C*log(h)/h)/(1-tau*c1*h)*fval
[~,fval] = fminbnd(@(x) -abs(splinef(x)),a,b)
eBoundEst = -tau*C*h^2*fval/(1-tau*C*h^2)



