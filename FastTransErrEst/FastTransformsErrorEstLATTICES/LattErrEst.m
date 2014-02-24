%Automatic cubature with Sobol sequences

%% Garbage collection and initialization
format compact %remove blank lines from output
format long e %lots of digits
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',30,'defaulttextfontsize',30) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %latex axis labels
tic

%% Initialize parameters
mmax=20; %maximum number of points is 2^mmax
mmin=6; %initial number of points is 2^mmin
mlag=5;
%testfun=@(x) x; exactinteg=1/2; d=1; %test function
%testfun=@(x) x.^2; exactinteg=1/3; d=1; %test function
a=20; testfun=@(x) sin(a*x); exactinteg=(1-cos(a))/a; d=1; %test function
%testfun=@(x) x(:,1).*x(:,2); exactinteg=1/4; d=2; %test function
%testfun=@(x) sin(x(:,1)).*x(:,2)+exp(x(:,1)); exactinteg=(1-cos(1))/2 + (exp(1)-1); d=2; %test function
sobstr=sobolset(d);
sobstr=scramble(sobstr,'MatousekAffineOwen');
sobol=qrandstream(sobstr);
Stilde=zeros(mmax-mmin+1,1);
appxinteg=zeros(mmax-mmin+1,1);


%% Initial points and FWT
n=2^mmin;
xpts=rand(sobol,n,d);
y=testfun(xpts);
%keyboard
yval=y;
yfwt=fwht(y);

%% Compute initial FWT
for l=0:mmin-1
   nl=2^l;
   nmminlm1=2^(mmin-l-1);
   ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+oddval)/2;
   y(~ptind)=(evenval-oddval)/2;
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
nllstart=2^(mmin-mlag-1);
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));

%disp('line 68')

%% Loop over m
for m=mmin+1:mmax 
   %disp('line 72')
   n=2^m;
   mnext=m-1;
   nnext=2^mnext;
   xnext=rand(sobol,nnext,d);
   ynext=testfun(xnext);
   yval=[yval; ynext];

   %% Compute initial FWT on next points
   for l=0:mnext-1
      nl=2^l;
      nmminlm1=2^(mnext-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      evenval=ynext(ptind);
      oddval=ynext(~ptind);
      ynext(ptind)=(evenval+oddval)/2;
      ynext(~ptind)=(evenval-oddval)/2;
   end
   %disp('line 90')

   %% Compute FWT on all points
   y=[y;ynext];
   nl=2^mnext;
   ptind=[true(nl,1); false(nl,1)];
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+oddval)/2;
   y(~ptind)=(evenval-oddval)/2;
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
   nllstart=2^(m-mlag-1);
   Stilde(m-mmin+1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));

   %% Approximate integral
   appxinteg(m-mmin+1)=mean(yval);

end

%% Error 
trueerr=abs(exactinteg-appxinteg);
%disp([appxinteg trueerr Stilde])

%% Plot results
figure
minexp=floor(mmin*log10(2));
maxexp=ceil(mmax*log10(2));
h=loglog(2.^(mmin:mmax),trueerr,'b.',2.^(mmin:mmax),Stilde,'rs');
set(h(2),'MarkerFaceColor','r','MarkerSize',10);
set(gca,'Xtick',10.^(minexp:maxexp))
axis([10.^[minexp maxexp],1e-13 1])

toc

