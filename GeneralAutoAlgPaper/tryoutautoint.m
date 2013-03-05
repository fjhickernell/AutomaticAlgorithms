%try out adapttrap

%% Garbage cleanup
format long e
clear all
close all
clc

%% Program parameters
param.tau=1000; %size of the cone
param.nmin=1000; %minimum number of trapezoids
param.nmax=1e6; %maximum number of trapezoids
param.tol=1e-8; %error tolerance

%% Simulation parameters
nrep=10000;

avec=exp((3*rand(nrep,1)+1).*log(10));
% cvec=10.^(1+5*(rand(nrep,1).^5));   %a

muvec=1./(sqrt(2).*avec)+rand(nrep,1).*(1-sqrt(2)./avec);
% zvec=rand(nrep,1);   %mu
tauvec=[10 100 1000]; %tau
nM=3;

%% Simulation
ratiovec=zeros(nrep,1);
Qmat=zeros(nrep,nM);
ntrapmat=Qmat;
errestmat=Qmat;
timemat=Qmat;
exactintmat=Qmat;

for i=1:nrep
    %% Integrand
    a=avec(i);
    mu=muvec(i);
    
%     f=@(x) exp(-(c*(x-z)).^2).*(2*c./...
%         (sqrt(pi)*(erf(c*(1-z))+erf(c*z))));    
    f=@(x) exp(-(a*(x-mu)).^2);

    exactint=0.5*sqrt(pi)*(erf(a*mu)-erf(a*(mu-1)))/a;
    exactintmat(i,:)=exactint;

%     param.exactinteg=1;
    
%     fp1norm=2 - exp(-(c*(z-1)).^2) - exp(-(c*z)^2);
%     fpp1norm=2*c*(2*sqrt(2/exp(1)) - c*((1-z)*exp(-(c*(z-1)).^2)...
%         -z*exp(-(c*z)^2)));

    fp1norm=-(exp(-(a.*(mu-1)).^2) + exp(-(a.*mu).^2)-2);
    fpp1norm=2*a.*(2*sqrt(2/exp(1))-a.*((1-mu).*exp(-...
        (a.*(mu-1)).^2)- a.*mu.*exp(-(a.*mu).^2)));

%     fpsupnornm=c*sqrt(2/exp(1));
%     fppsupnorm=2*c^2;

    param.ratio=fpp1norm/fp1norm;
%     param.ratio=fppsupnorm/fpsupnorm;
    ratiovec(i)=param.ratio;

    %%Integrate
    for j=1:nM
        param.tau=tauvec(j);
        param.nmin=1;
        tic
        [Q,param]=autoint(f,param);
        timemat(i,j)=toc;
        Qmat(i,j)=Q;
        ntrapmat(i,j)=param.ntrap;
        errestmat(i,j)=param.errest;
    end
end

trueerrormat=abs(exactintmat-Qmat);
successrate=mean(trueerrormat<=param.tol,1);
ratiomat=repmat(ratiovec,1,nM);
taumat=repmat(tauvec,nrep,1);
ratiosmaller=mean(ratiomat<=taumat,1);
violaterate=mean((trueerrormat>param.tol)&(ratiomat<=taumat),1);
avgtime=mean(timemat,1);

save('TrapezoidGauss.mat')

display(successrate);