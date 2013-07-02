%try out adapttrap

%% Garbage cleanup
format long e
tstart=tic;

%% Program parameters
param.nmax=1e8; %maximum number of trapezoids
param.tol=1e-8; %error tolerance
paramnon.ntrap=1e6;

%% Simulation parameters
nrep=10000;
%cvec=10.^(1+5*(rand(nrep,1).^5));
cvec=10.^(1+3*(rand(nrep,1)));
zvec=1./(sqrt(2)*cvec)+ (1-sqrt(2)./cvec).*rand(nrep,1);
Mvec=[10 100 1000];
nM=length(Mvec);

%% Simulation
ratiovec=zeros(nrep,1);
Qmat=zeros(nrep,nM);
ntrapmat=Qmat;
errestmat=Qmat;
timemat=Qmat;
newtaumat=Qmat;
Qquadvec=ratiovec;
Qquadgkvec=ratiovec;
Qchebvec=ratiovec;
Qnonadaptvec=ratiovec;
quadtimevec=ratiovec;
quadgktimevec=ratiovec;
chebtimevec=ratiovec;
nonadapttimevec=ratiovec;
for i=1:nrep
    if floor(i/50)==i/50, disp(i), end
    
    %% Integrand
    c=cvec(i);
    z=zvec(i);
    f=@(x) exp(-(c*(x-z)).^2).*(2*c./...
        (sqrt(pi)*(erf(c*(1-z))+erf(c*z))));
    param.exactinteg=1;
    fp1norm=2 - exp(-(c*(z-1)).^2) - exp(-(c*z)^2);
    fpp1norm=2*c*(2*sqrt(2/exp(1)) - c *((1 - z)*exp(-(c*(z-1)).^2)...
        - z*exp(-(c*z)^2)));
    param.ratio=fpp1norm/fp1norm;
    ratiovec(i)=param.ratio;

    %%Integrate
    for j=1:nM
        param.tau=Mvec(j);
        param.nmin=1;
        tic
        [Q,param]=adapttrapcheck(f,param);
        timemat(i,j)=toc;
        Qmat(i,j)=Q;
        ntrapmat(i,j)=param.ntrap;
        errestmat(i,j)=param.errest;
        newtaumat(i,j)=param.tau;
    end
    
    tic;
    Qquadvec(i)=quad(f,0,1,param.tol);
    quadtimevec(i)=toc;
    tic;
    Qquadgkvec(i)=quadgk(f,0,1,'AbsTol',param.tol);
    quadgktimevec(i)=toc;
    tic;
    Qchebvec(i)=sum(chebfun(f,[0,1]));
    chebtimevec(i)=toc;
    tic;
    Qnonadaptvec(i)=nonadapttrap(f,paramnon);
    nonadapttimevec(i)=toc;
end

%% Postprocessing
trueerrormat=abs(1-Qmat);
successrate=mean(trueerrormat<=param.tol,1)
ratiomat=repmat(ratiovec,1,nM);
Mmat=repmat(Mvec,nrep,1);
ratiosmaller=mean(ratiomat<=Mmat,1)
rationewsmaller=mean(ratiomat<=newtaumat,1)
violaterate=mean((trueerrormat>param.tol)&(ratiomat<=newtaumat),1)
%avgtime=mean(timemat,1)
avgtime=geo_mean(timemat,1)
quadtrueerrorvec=abs(1-Qquadvec);
quadsuccessrate=mean(quadtrueerrorvec<=param.tol,1)
%quadavgtime=mean(quadtimevec,1)
quadavgtime=geo_mean(quadtimevec,1)
quadgktrueerrorvec=abs(1-Qquadgkvec);
quadgksuccessrate=mean(quadgktrueerrorvec<=param.tol,1)
%quadgkavgtime=mean(quadgktimevec,1)
quadgkavgtime=geo_mean(quadgktimevec,1)
chebtrueerrorvec=abs(1-Qchebvec);
chebsuccessrate=mean(chebtrueerrorvec<=param.tol,1)
%chebavgtime=mean(chebtimevec,1)
chebavgtime=geo_mean(chebtimevec,1)
nonadapttrueerrorvec=abs(1-Qnonadaptvec);
nonadaptsuccessrate=mean(nonadapttrueerrorvec<=param.tol,1)
nonadaptavgtime=mean(nonadapttimevec,1)

%% Saving data
save(['TrapezoidGauss' datestr(now,'dd-mmm-yyyy-HH:MM:SS') '.mat'])
toc(tstart)
