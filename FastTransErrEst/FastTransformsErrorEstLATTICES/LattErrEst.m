%Automatic cubature with lattices

%% Garbage collection and initialization
format compact %remove blank lines from output
format long e %lots of digits
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',15,'defaulttextfontsize',15) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',30) %latex axis labels
tic

%% Particular example where the Fourier cone is bounded by b^-k (check paper on function approx)
f=@(b,x) (b^2-1)./(b^2+1-2*b*cos(2*pi()*x));
g=@(x) 8./(10-6*cos(2*pi()*x));
%% Bernoulli polynomial
ber=@(x) x.^2-x+1/4;

%% Initialize parameters
mmax=15; %maximum number of points is 2^mmax
mmin=6; %initial number of points is 2^mmin
mlag=5;
latticeseq_b2('init0'); %initializing lattice numbers generator
%testfun=@(x) x; exactinteg=1/2; d=1; %test function
%testfun=@(x) x.^2; exactinteg=1/3; d=1; %test function
%testfun=@(x) g(x); exactinteg=1; d=1; %test function
%testfun=@(x) ber(x); exactinteg=1/12; d=1; %test function
%testfun=@(x) exp(cos(2*pi*x)); exactinteg=1.266065877752008; d=1; %test function
%a=20; testfun=@(x) sin(a*x); exactinteg=(1-cos(a))/a; d=1; %test function
%testfun=@(x) 2*x.*(x<0.5)+(2-2*x).*(0.5<=x); exactinteg=1/2; d=1; %test function
%testfun=@(x) (x-0).^2.*(x-1).^2; exactinteg=1/30; d=1; %test function
%testfun=@(x) (x(:,1)-0).^2.*(x(:,1)-1).^2.*(x(:,2)-0).^2.*(x(:,2)-1).^2; exactinteg=1/30^2; d=2; %test function
%testfun=@(x) x(:,1).*x(:,2); exactinteg=1/4; d=2; %test function
%testfun=@(x) g(x(:,1)).*g(x(:,2)); exactinteg=1; d=2; %test function
%testfun=@(x) ber(x(:,1)).*ber(x(:,2)); exactinteg=1/12^2; d=2; %test function
testfun=@(x) sin(x(:,1)).*x(:,2)+exp(x(:,1)); exactinteg=(1-cos(1))/2 + (exp(1)-1); d=2; %test function
Stilde=zeros(mmax-mmin+1,1);
appxinteg=zeros(mmax-mmin+1,1);

%% Adding the Baker's transform? For 1D
%testfun=@(x) testfun(1-2*abs(x-1/2)); exactinteg=exactinteg; d=1; %test function with Baker transform
%% Adding the Baker's transform? For 2D
%testfun=@(x) testfun([1-2*abs(x(:,1)-1/2),1-2*abs(x(:,2)-1/2)]); exactinteg=exactinteg; d=2; %test function with Baker transform

%% Initial points and FFT
n=2^mmin;
xpts=latticeseq_b2(d,n)';
y=testfun(xpts);
%keyboard
yval=y;
%yfft=fft(y);

%% Compute initial FFT
for l=0:mmin-1
   nl=2^l;
   nmminlm1=2^(mmin-l-1);
   ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
   coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(2*nl))';
   coefv=repmat(coef,nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+coefv.*oddval)/2;
   y(~ptind)=(evenval-coefv.*oddval)/2;
end

%% Approximate integral
appxinteg(1)=mean(yval);

%% Create kappanumap
kappanumap=(1:n)'; %initialize map
for l=mmin-1:-1:1
   nl=2^l;
   oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
   newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
   flip=find(newone>oldone); %
   temp=kappanumap(nl+1+flip);
   kappanumap(nl+1+flip)=kappanumap(1+flip);
   kappanumap(1+flip)=temp;
   %keyboard
end
%[y y(kappanumap)];
%keyboard

%% Compute Stilde
nllstart=int64(2^(mmin-mlag-1)); %int64 to convert to integer
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));

%disp('line 68')

%% Loop over m
for m=mmin+1:mmax 
   %disp('line 72')
   n=2^m;
   mnext=m-1;
   nnext=2^mnext;
   xnext=latticeseq_b2(d,nnext)';
   ynext=testfun(xnext);
   yval=[yval; ynext];

   %% Compute initial FFT on next points
   for l=0:mnext-1
      nl=2^l;
      nmminlm1=2^(mnext-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(2*nl))';
      coefv=repmat(coef,nmminlm1,1);
      evenval=ynext(ptind);
      oddval=ynext(~ptind);
      ynext(ptind)=(evenval+coefv.*oddval)/2;
      ynext(~ptind)=(evenval-coefv.*oddval)/2;
end

   %disp('line 90')

   %% Compute FFT on all points
   y=[y;ynext];
   nl=2^mnext;
   ptind=[true(nl,1); false(nl,1)];
   coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(2*nl))';
   coefv=repmat(coef,nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+coefv.*oddval)/2;
   y(~ptind)=(evenval-coefv.*oddval)/2;

   %disp('line 100')

   %% Update kappanumap
   kappanumap=[kappanumap; (nnext+1:n)']; %initialize map
   for l=m-1:-1:1
      nl=2^l;
      oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
      newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
      flip=find(newone>oldone); %
      temp=kappanumap(nl+1+flip);
      kappanumap(nl+1+flip)=kappanumap(1+flip);
      kappanumap(1+flip)=temp;
   end
   %disp('line 114')

   %[y y(kappanumap)];
   %keyboard

   %% Compute Stilde
   nllstart=int64(2^(m-mlag-1));
   Stilde(m-mmin+1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   
   
   %% Approximate integral
   appxinteg(m-mmin+1)=mean(yval);

end
figure
ymap=y(kappanumap);
loglog(abs(ymap))


%% Error 
trueerr=abs(exactinteg-appxinteg);
%disp([appxinteg trueerr Stilde])
errorstilde=Stilde.*2.^(-mmin:-1:-mmax)';

%% Plot results
figure
minexp=floor(mmin*log10(2));
maxexp=ceil(mmax*log10(2));
% h=loglog(2.^(mmin:mmax),trueerr,'b.',2.^(mmin:mmax),Stilde,'rs');
% set(h(2),'MarkerFaceColor','r','MarkerSize',10);
% set(gca,'Xtick',10.^(minexp:maxexp))
% axis([10.^[minexp maxexp],1e-17 10])
figure
h=loglog(2.^(mmin:mmax),trueerr,'b.',2.^(mmin:mmax),errorstilde,'rs');
set(h(2),'MarkerFaceColor','r','MarkerSize',10);

toc