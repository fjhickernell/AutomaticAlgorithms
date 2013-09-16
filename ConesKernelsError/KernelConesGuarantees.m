%Cones and Kernel Approximation for Guarantees
%% Preliminaries
clear all
close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
testfun=@(x) sin(8*pi*x);
tau=8*pi;
epsilon = 1e-2;

%% Kernel
a=10;
kernel=@(x,t) exp(-(a.^2)*(bsxfun(@minus,x,t')).^2);

%% Data and spline approximation
n=25;
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
error=sqrt(mean((testfun(xtest)-splinef(xtest)).^2))

normHsplinef=sqrt(c'*y)
xpxt=bsxfun(@plus,xnode,xnode');
xmxt=bsxfun(@minus,xnode,xnode');
Ktildemat=(sqrt(pi/2)/(2*a)) ...
    *exp(-((a^2)/2)*(xmxt.^2)) ...
    .*(-erf((a/sqrt(2))*(xpxt-2)) + erf((a/sqrt(2))*xpxt));
condKtilde=cond(Ktildemat)
Herrbd=sqrt(1-trace(Ktildemat/Kmat))
guesserrest=Herrbd*normHsplinef

%% Algorithm 1 Stage 1
xtxt=bsxfun(@times,xnode,xnode');
xpxtsq=bsxfun(@plus,xnode.^2,xnode'.^2);
Htildemat=1/4*a*(exp(-a^2*(2+xpxtsq+xtxt-2*xpxt)) ...
    .*(2*a*exp(a^2*xtxt).*(-2+xpxt)+exp(1/2*a^2*(bsxfun(@plus,xnode.^2-4*xnode,(xnode'-2).^2)+4*xtxt)) ...
    *sqrt(2*pi).*(-1+a^2*xmxt.^2)*erf(a*(-2+xpxt)/sqrt(2))) ...
    -exp(-a^2*(xpxtsq+xtxt)).*(2*a*exp(a^2*xtxt).*xpxt+exp(1/2*a^2*(xpxtsq+4*xtxt)) ...
    *sqrt(2*pi)*(-1+a^2*xmxt.^2).*erf(a*xpxt/sqrt(2))));
condHtilde=cond(Htildemat)
normDsplinef=sqrt(c'*Htildemat*c);
htilde=sqrt(2*a^2-trace(Htildemat/Kmat))