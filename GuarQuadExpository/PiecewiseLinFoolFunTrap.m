% Piecewise liear fooling function for the trapezoid rule

%% Garbage collection
format compact %remove blank lines from output
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'DefaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels

%% Parameters
n=12; N=2^n; %number of nodes - 1 to define function
m=6; M=2^m; %number of trapezoids in finer rule
mmin1=m-1; Mov2=2^mmin1; %number of traepzoids in coarser rule
L=2^(n-m);
twoL=2*L;
f0=0;

%% Compute s and d values
jnodevec1=1:N-1;
integsvec=(N-jnodevec1).*(N-1-jnodevec1)/(2*N);
integsvec=integsvec ...
    -(N-1:-1:1)*((N-1)/(2*N));

jtfinevec1=1:N-L;
cljtfinevec=ceil(jtfinevec1/L);
tfinesvec=(M-cljtfinevec).*...
    (N-L+L*cljtfinevec-2*jtfinevec1)/(2*M);
tfinesvec=[tfinesvec zeros(1,L-1)];
tfinesvec=tfinesvec ...
    -(N-1:-1:1)*((N-L)/(2*N));

jtcoarsevec1=1:N-twoL;
cljtcoarsevec=ceil(jtcoarsevec1/twoL);
tcoarsesvec=(Mov2-cljtcoarsevec).*...
    (N-twoL+twoL*cljtcoarsevec-2*jtcoarsevec1)/M;
tcoarsesvec=[tcoarsesvec zeros(1,twoL-1)];
tcoarsesvec=tcoarsesvec ...
    -(N-1:-1:1)*((N-twoL)/(2*N));

Bmat=[integsvec-tfinesvec; ... %error
    tfinesvec-tcoarsesvec; ... %error est
    integsvec]'; %integral
%    N-1:-1:1]'; %periodicity

% Cmat=[(L-1)/2; ... %error
%     L/2; ... %error estimate
%     N]'; %periodicity

dvec=[1; ... %error
    3*0; ... %error estimate
    1]; %integral
%    0]; %periodicity

[U2,Lamb2,V2]=svd(Bmat,0);
%[U3,Lamb3,V3]=svd(Lamb2\(V2'*Cmat'),0);
pc1vec=Lamb2\(V2'*dvec);
%pc2vec=U3'*pc1vec;
%yvec=V3*(Lamb3\pc2vec);
%s=U2*(pc1vec-U3*pc2vec);
s=U2*pc1vec;
d0=-((N-1:-1:1)*s)/N;
%d=cumsum([yvec; s]);
d=cumsum([d0; s]);
f=cumsum([f0; d]);

%Plot resuts so far
xfun=(0:N)'/N;
xftrap=(0:L:N)'/N;
xctrap=(0:twoL:N)'/N;
fftrap=f(1:L:N+1);
fctrap=f(1:twoL:N+1);
plot(xfun,f,'k-',...
    xftrap,fftrap,'b-', ...
    xctrap,fctrap,'r-')

figure;
plot(xfun(1:2*twoL),f(1:2*twoL),'k-',...
    xftrap(1:5),fftrap(1:5),'b-', ...
    xctrap(1:3),fctrap(1:3),'r-')


%Check results
integ=mean(f(1:N));
finetrap=mean(fftrap(1:M));
coarsetrap=mean(fctrap(1:Mov2));
err=abs(integ-finetrap);
errest=abs(finetrap-coarsetrap)/3;
disp(['      Integral = ' num2str(integ)])
disp(['         Error = ' num2str(err)])
disp(['Error Estimate = ' num2str(errest)]);
    


