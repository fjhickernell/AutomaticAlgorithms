%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
testfun=@(x) sin(3*pi*x);
tau=1000;
epsilon = 1e-2;

%% Kernel
a=2;
kernel=@(x,t) exp(-(a.^2)*(bsxfun(@minus,x,t')).^2);

%% Data and spline approximation
n = 5
h = 1/n;
xnode=linspace(0,1,n)';
Kmat=kernel(xnode,xnode);
condK=cond(Kmat)
y=testfun(xnode);
c=Kmat\y;
splinef=@(x) kernel(x,xnode)*c;

%% Plot of test function and spline
ntest=5000;
xtest=linspace(0,1,ntest)';
ytest = interp1q(xnode,y,xtest);
plot(xtest,ytest,'b-',...
    xtest,splinef(xtest),'r--',...
    xnode,y,'k.','markersize',20,...
    'linewidth',3)

%% Root mean square error
error2=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))
error_sup=max(abs(testfun(xtest)-splinef(xtest)))
% normHsplinef=sqrt(c'*y)
% xpxt=bsxfun(@plus,xnode,xnode');
% xmxt=bsxfun(@minus,xnode,xnode');
% Ktildemat=(sqrt(pi/2)/(2*a)) ...
%     *exp(-((a^2)/2)*(xmxt.^2)) ...
%     .*(-erf((a/sqrt(2))*(xpxt-2)) + erf((a/sqrt(2))*xpxt));
% condKtilde=cond(Ktildemat)
% Herrbd=sqrt(1-trace(Ktildemat/Kmat))
% guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
Ktestmat = kernel(xtest,xnode);
normDsplinef = norm(Ktestmat*c,Inf);
GaussConst = 1;
eBoundEst = tau*sqrt(exp(-GaussConst/h))/.9*normDsplinef

% xtxt=bsxfun(@times,xnode,xnode');
% xpxtsq=bsxfun(@plus,xnode.^2,xnode'.^2);
% Htildemat=1/4*a*(exp(-a^2*(2+xpxtsq+xtxt-2*xpxt)) ...
%     .*(2*a*exp(a^2*xtxt).*(-2+xpxt)+exp(1/2*a^2*(bsxfun(@plus,xnode.^2-4*xnode,(xnode'-2).^2)+4*xtxt)) ...
%     *sqrt(2*pi).*(-1+a^2*xmxt.^2)*erf(a*(-2+xpxt)/sqrt(2))) ...
%     -exp(-a^2*(xpxtsq+xtxt)).*(2*a*exp(a^2*xtxt).*xpxt+exp(1/2*a^2*(xpxtsq+4*xtxt)) ...
%     *sqrt(2*pi)*(-1+a^2*xmxt.^2).*erf(a*xpxt/sqrt(2))));
% condHtilde=cond(Htildemat)
% normDsplinef=sqrt(c'*Htildemat*c);
% htilde=sqrt(2*a^2-trace(Htildemat/Kmat))