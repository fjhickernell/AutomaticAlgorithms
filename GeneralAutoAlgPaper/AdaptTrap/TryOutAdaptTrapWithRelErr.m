%try out adapttrap

%% Garbage cleanup
format long e

%% Program parameters
param.Mcone=1000; %size of the cone
param.nmin=1; %minimum number of trapezoids
param.nmax=1e8; %maximum number of trapezoids
%param.tol=1e-8; %error tolerance
param.abstol=1e-6; %absolute error tolerance
param.reltol=1e-6; %relative error tolerance
param.theta=0;
totalerr=(1-param.theta)*param.abstol+param.theta*param.reltol;
paramnon.ntrap=1e6;

%% Simulation parameters
nrep=1000; %Number of replications
cvec=10.^(1+5*(rand(nrep,1).^5));
zvec=rand(nrep,1);
exactinteg=10.^(-1+2*rand(nrep,1));
Mvec=[10 100 1000];
nM=3;

%% Simulation
ratiovec=zeros(nrep,1);
Qmat=zeros(nrep,nM);
ntrapmat=Qmat;
errestmat=Qmat;
timemat=Qmat;
Qquadvec=ratiovec;
Qquadgkvec=ratiovec;
Qchebvec=ratiovec;
Qnonadaptvec=ratiovec;
quadtimevec=ratiovec;
quadgktimevec=ratiovec;
chebtimevec=ratiovec;
nonadapttimevec=ratiovec;
for i=1:nrep
    if  floor(i/10)==i/10, disp(i), end
    %% Integrand
    c=cvec(i);
    z=zvec(i);
    f=@(x) exactinteg(i)*exp(-(c*(x-z)).^2).*(2*c./...
        (sqrt(pi)*(erf(c*(1-z))+erf(c*z))));
    param.exactinteg=1;
    fp1norm=2 - exp(-(c*(z-1)).^2) - exp(-(c*z)^2);
    fpp1norm=2*c*(2*sqrt(2/exp(1)) - c *((1 - z)*exp(-(c*(z-1)).^2) - z*exp(-(c*z)^2)));
    param.ratio=fpp1norm/fp1norm;
    ratiovec(i)=param.ratio;

    %%Integrate
    for j=1:nM
        param.Mcone=Mvec(j); %parameter defining the cone
        param.nmin=1; %minimum number of trapezoids
        tic
        [Q,param]=adapttrapwrel(f,param);
        timemat(i,j)=toc;
        Qmat(i,j)=Q;
        ntrapmat(i,j)=param.ntrap;
        errestmat(i,j)=param.errest;
    end
    
    if param.theta==0;
        tic;
        Qquadvec(i)=quad(f,0,1,param.abstol);
        quadtimevec(i)=toc;    
    end
    if param.theta==0;
        tic;
        Qquadgkvec(i)=quadgk(f,0,1,'AbsTol',param.abstol,'RelTol',0);
        quadgktimevec(i)=toc;
    else
        tic;
        Qquadgkvec(i)=quadgk(f,0,1,'AbsTol',0,'RelTol',param.reltol);
        quadgktimevec(i)=toc;
    end
    tic;
    Qchebvec(i)=sum(chebfun(f,[0,1]));
    chebtimevec(i)=toc;
    tic;
    Qnonadaptvec(i)=nonadapttrap(f,paramnon);
    nonadapttimevec(i)=toc;
end
trueerrormat=abs(repmat(exactinteg,1,nM)-Qmat);
successrate=mean(trueerrormat<=totalerr,1)
ratiomat=repmat(ratiovec,1,nM);
Mmat=repmat(Mvec,nrep,1);
ratiosmaller=mean(ratiomat<=Mmat,1)
violaterate=mean((trueerrormat>totalerr)&(ratiomat<=Mmat),1)
avgtime=mean(timemat,1)
quadtrueerrorvec=abs(exactinteg-Qquadvec);
quadsuccessrate=mean(quadtrueerrorvec<=totalerr,1)
quadavgtime=mean(quadtimevec,1)
quadgktrueerrorvec=abs(exactinteg-Qquadgkvec);
quadgksuccessrate=mean(quadgktrueerrorvec<=totalerr,1)
quadgkavgtime=mean(quadgktimevec,1)
chebtrueerrorvec=abs(exactinteg-Qchebvec);
chebsuccessrate=mean(chebtrueerrorvec<=totalerr,1)
chebavgtime=mean(chebtimevec,1)
nonadapttrueerrorvec=abs(exactinteg-Qnonadaptvec);
nonadaptsuccessrate=mean(nonadapttrueerrorvec<=totalerr,1)
nonadaptavgtime=mean(nonadapttimevec,1)
save(['TrapezoidGauss' datestr(now,'dd-mmm-yyyy-HH-MM-SS') '.mat'])

