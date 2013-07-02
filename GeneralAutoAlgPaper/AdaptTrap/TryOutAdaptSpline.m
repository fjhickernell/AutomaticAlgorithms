%try out adapttrap

%% Garbage cleanup
format long e

%% Program parameters
param.nmin=1000; %minimum number of intervals
param.nmax=1e8; %maximum number of intervals
param.tol=1e-8; %error tolerance

%% Simulation parameters
nrep=100;
%cvec=10.^(1+5*(rand(nrep,1).^5));
%cvec=10.^(1+3*(rand(nrep,1)));
cvec=10.^(4*(rand(nrep,1)));
%zvec=1./(sqrt(2)*cvec)+ (1-sqrt(2)./cvec).*rand(nrep,1);
zvec=rand(nrep,1);
Mvec=[10 100 1000];
nM=3;

%% Simulation
ratiovec=zeros(nrep,1);
nintmat=zeros(nrep,nM);
errestmat=nintmat;
errmat=nintmat;
timemat=nintmat;
ntest1=1000;
ntest2=1000;
for i=1:nrep
    if floor(i/10)==i/10, disp(i), end
    
    %% Function to be approximated
    c=cvec(i);
    z=zvec(i);
    f=@(x) exp(-(c*(x-z)).^2);
    fpinfnorm=c*sqrt(2/exp(1));
    fppinfnorm=2*c*c;
    param.ratio=fppinfnorm/fpinfnorm;
    ratiovec(i)=param.ratio;
    
    %% Sample points
    xtest1=z-erfinv(erf(c*z)+((1:ntest1-1)/ntest1)*(erf(c*(z-1))-erf(c*z)))/c;
    xtest2=(0:ntest2)/ntest2;
    xtest=sort([xtest1 xtest2]);
    ftest=f(xtest);

    %% Approximate
    for j=1:nM
        param.tau=Mvec(j);
        param.nmin=1;
        tic
        [fspline,param]=adaptlinspline(f,param);
        timemat(i,j)=toc;
        nintmat(i,j)=param.nint;
        errestmat(i,j)=param.errest;
        errmat(i,j)=max(abs(ftest-fspline(xtest)));
    end
    
end
successrate=mean(errmat<=param.tol,1)
ratiomat=repmat(ratiovec,1,nM);
Mmat=repmat(Mvec,nrep,1);
ratiosmaller=mean(ratiomat<=Mmat,1)
violaterate=mean((errmat>param.tol)&(ratiomat<=Mmat),1)
avgtime=mean(timemat,1)
save(['LinearSpline' datestr(now,'dd-mmm-yyyy-HH:MM:SS') '.mat'])

