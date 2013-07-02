%Try out adaptspline
close all, clear all
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)

%% Program parameters
param.tau=100; %size of the cone
param.nmin=1000; %minimum number of trapezoids
param.nmax=1e8; %maximum number of trapezoids
param.tol=1e-8; %error tolerance

testf=@(x) sin(10*x);
tic
[fspline,param]=adaptlinspline(testf,param)
toc
x=sort(rand(4000,1)); 
figure;
plot(x,testf(x),'b',x,fspline(x),'k','linewidth',2)
figure;
error=testf(x)-fspline(x);
finfnorm=max(error)
semilogy(x,abs(testf(x)-fspline(x)),'b','linewidth',2)