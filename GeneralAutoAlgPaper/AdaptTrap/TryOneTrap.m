%Try out adaptspline
close all, clear all
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)

%% Program parameters
param.tau=10; %size of the cone
oldtau=param.tau;
param.nmin=10; %minimum number of trapezoids
param.nmax=1e8; %maximum number of trapezoids
param.tol=1e-8; %error tolerance

%a=10; z=rand(1); testf=@(x) sin(a*(x-z)); exactinteg=(cos(a*z)-cos(a*(1-z)))/a;
a=1e1; z=rand(1); testf=@(x) exp(-(a*(x-z)).^2).*(2*a./(sqrt(pi)*(erf(a*(1-z))+erf(a*z)))); exactinteg=1;
%a=100; testf=@(x) a*x; exactinteg=a/2;
testf
tic
[Q,param]=adapttrapcheck(testf,param);
toc
err=exactinteg-Q;
disp(['Exact integral = ' num2str(exactinteg)])
disp(['Approximate integral = ' num2str(Q)])
disp(['Error = ' num2str(err)])
disp(['Computational Cost = ' int2str(param.ntrap+1)])
if param.success; succstr='true'; else succstr='false'; end
disp(['Success = ' succstr])
if param.tau == oldtau
    disp('No change in tau')
else
    disp(['tau changed from ' num2str(oldtau) ' to ' num2str(param.tau)])
end
