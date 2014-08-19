%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
n = 500
a = 0; b = 1;
xnode=linspace(a,b,n)';
h = (b-a)/(2*n); R = (b-a)/2;
testfun=@(x) sin(.1*x);
gamma = 1;
kernel=@(x,t) exp(-gamma*(bsxfun(@minus,x,t')).^2);
tau = 2;
c0 = R/6; C1 = 2304*gamma; C2 = log(4)+16; C3 = C1*exp(C2);
%C = min(c0/2,1/(exp(1)*C3))/2;
C = min(c0/2,1/sqrt(exp(1)*C3))/2;
%C4 = 256/3; C5 = 1/2; c1 = 2*gamma^2*C5*C4^3*(6+C5*C4);
%h0 = min([1/(tau*c1) c0 R/C4 sqrt(3)/C4/sqrt(2*gamma)]);
if tau<=1
    h0 = c0
else
    h0 = min([c0 C/(2*log(tau))])
end
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
error2=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))
error_sup=max(abs(testfun(xtest)-splinef(xtest)))
normHsplinef=sqrt(c'*y)
% xpxt=bsxfun(@plus,xnode,xnode');
% xmxt=bsxfun(@minus,xnode,xnode');
% Ktildemat=(sqrt(pi/2)/(2*a)) ...
%     *exp(-((a^2)/2)*(xmxt.^2)) ...
%     .*(-erf((a/sqrt(2))*(xpxt-2)) + erf((a/sqrt(2))*xpxt));
% condKtilde=cond(Ktildemat)
% Herrbd=sqrt(1-trace(Ktildemat/Kmat))
% guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
%[~,fval] = fminbnd(@(x) -FcnToMax(x,xnode,c,gamma,1),a,b);
%eBoundEst = -tau*exp(C*log(h)/h)/(1-tau*c1*h)*fval
%[~,fval] = fminbnd(@(x) -abs(splinef(x)),a,b);
%eBoundEst = -tau*exp(C*log(h)/h)/(1-tau*exp(C*log(h)/h))*fval
[~,fval] = fminbnd(@(x) -abs(splinef(x)),a,b)
eBoundEst = -tau*exp(-C/(2*h))/(1-tau*exp(-C/(2*h)))*fval


