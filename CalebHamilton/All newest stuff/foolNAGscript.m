%This script is being set up to specifically fool the NAG toolbox functions
clear all
format long
format compact
naginfo.filename='scriptValues.mat'; %This file name will be passed throughout the script
xsample=[];
save(naginfo.filename,'xsample')
%%
%------------------------Inputs-----------------------------
naginfo.morepeaks=false;
oneplot=true;
naginfo.lower=0;
naginfo.upper=1;
naginfo.c=10;
naginfo.p=5;
naginfo.coefficient=100;
naginfo.degree=1;
naginfo.sign=1;
naginfo.RegFunc=@(x) naginfo.coefficient.*x.^naginfo.degree;
naginfo.RegFuncprime=... %and its derivative
    @(x) naginfo.degree.*naginfo.coefficient*x.^(naginfo.degree-1); %slowly varying
naginfo.RegFuncdubprime=... %and its second derivative
    @(x) naginfo.degree.*(naginfo.degree-1).*naginfo.coefficient*x.^(naginfo.degree-2);
naginfo.quadtype='d01ah';
%%
peaks=@(x) peakyfunction(x,naginfo);
save('nagPeaks.mat','peaks','naginfo')
fname='snoopNAG';
switch naginfo.quadtype
    case 'd01ah'
        callautoalg=@(name,lower,upper) d01ah(lower,upper,1e-6,name,int32(100));
    case'd01aj'
        callautoalg=@(name,lower,upper) d01aj(name,lower,upper,1e-6,1e-6);
    case'd01ak'
        callautoalg=@(name,lower,upper) d01ak(name,lower,upper,1e-6,1e-6);
end
naginfo.original=callautoalg(fname,naginfo.lower,naginfo.upper);
naginfo.secondattempt=callautoalg('nagfoolpeaky',naginfo.lower,naginfo.upper);
%%
%--------------------------------Plots-------------------------------
xx=naginfo.lower:.005:naginfo.upper; %Evenly spaced points for plotting
[yplot,yprimeplot,ydubplot,naginfo]=peakyfunction(xx,naginfo);
if(~oneplot)
    figure(1)
    subplot(2,3,1),plot(xx,yplot,'b-',naginfo.sortedX,naginfo.RegFunc(naginfo.sortedX),'ro')
    subplot(2,3,3),plot(xx,yprimeplot),title('First Derivative')
    subplot(2,3,5),plot(xx,ydubplot),title('Second Derivative')
else
    figure(1)
    plot(xx,yplot,'b-',naginfo.sortedX,naginfo.RegFunc(naginfo.sortedX),'ro')
end
%%
%-----------------------------Integrals---------------------------
upperbnd=naginfo.sortedX(2:end);
lowerbnd=naginfo.sortedX(1:end-1);
withbumps=(upperbnd-lowerbnd)/2*(sqrt(pi)*gamma(naginfo.p+1)/gamma(1.5+naginfo.p));
failintegral=naginfo.secondattempt;
naginfo.realintegral=naginfo.c*(sum(withbumps))+naginfo.coefficient/(naginfo.degree+1)*...
    (naginfo.upper^(naginfo.degree+1)-naginfo.lower^(naginfo.degree+1));
%relative error
naginfo.error=abs((naginfo.realintegral-failintegral)/naginfo.realintegral);
%%
%--------------------------------Derivatives------------------------
%===========================Infinity Norm==============================
naginfo.secondzeros=[naginfo.sortedX(2:end)'; naginfo.sortedX(1:end-1)'; (.5*((naginfo.sortedX(2:end)...
    -naginfo.sortedX(1:end-1))/sqrt(2*naginfo.p-1)+naginfo.sortedX(1:end-1)+...
    naginfo.sortedX(2:end)))'; (.5*((-naginfo.sortedX(2:end)+naginfo.sortedX(1:end-1))/...
    sqrt(2*naginfo.p-1)+naginfo.sortedX(1:end-1)+naginfo.sortedX(2:end)))'];
[~,yp]=peaks(naginfo.secondzeros);
primemaxvec=max(abs(yp));
naginfo.thirdzeros=[naginfo.sortedX(2:end)'; naginfo.sortedX(1:end-1)'; ((naginfo.sortedX(1:end-1)+...
    naginfo.sortedX(2:end))/2)'; (.5*(sqrt(3)*(naginfo.sortedX(2:end)-naginfo.sortedX(1:end-1))/...
    sqrt(2*naginfo.p-1)+naginfo.sortedX(1:end-1)+naginfo.sortedX(2:end)))'; (.5*(sqrt(3)*(-naginfo.sortedX(2:end)...
    +naginfo.sortedX(1:end-1))/sqrt(2*naginfo.p-1)+naginfo.sortedX(1:end-1)+naginfo.sortedX(2:end)))'];
[~,~,ydub]=peaks(naginfo.thirdzeros);
dubmaxvec=max(abs(ydub));
primemax=max(primemaxvec);
dubmax=max(dubmaxvec);
naginfo.inormratio=dubmax/(primemax+naginfo.coefficient);

%==========================Two Norm========================
if(naginfo.morepeaks)
    prime2norm=max(1./sqrt(upperbnd-lowerbnd)*2*naginfo.p*pi^.25*sqrt(gamma(2*naginfo.p-1)/...
      gamma(2*naginfo.p+.5)));
    dub2norm=max(1./sqrt(upperbnd-lowerbnd).^3*2^(5/2)*naginfo.p*sqrt((3*sqrt(pi)*(naginfo.p-1)^2*...
     gamma(naginfo.p*2-3))/(gamma(2*naginfo.p-.5))));
    naginfo.ratio2norm=dub2norm/prime2norm;
else
    u=naginfo.sortedX(naginfo.largeststart(1)+1);
    l=naginfo.sortedX(naginfo.largeststart(1));
    prime2norm=max(1/sqrt(u-l)*2*naginfo.p*pi^.25*sqrt(gamma(2*naginfo.p-1)/...
    gamma(2*naginfo.p+.5)));
    dub2norm=max(1/sqrt(u-l).^3*2^(5/2)*naginfo.p*sqrt((3*sqrt(pi)*(naginfo.p-1)^2*...
      gamma(naginfo.p*2-3))/(gamma(2*naginfo.p-.5))));
    naginfo.ratio2norm=dub2norm/prime2norm;
end
%%
%----------------------------------Ratio(E)---------------------
E=.001:.005:1;
ratiofunc=64./((naginfo.sortedX(naginfo.largeststart+1)-naginfo.sortedX(naginfo.largeststart))^2*...
    (9*sqrt(3))*(3/2*sqrt(pi)*gamma(naginfo.p+1)/gamma(naginfo.p+3/2)*(1./E-1)+1));
figure(2)
plot(E,ratiofunc),title('R(E)'),xlabel('Relative Error'),ylabel('Infinity Norm Ratio')

%%
%----------------------------------Output---------------------
load(naginfo.filename)
naginfo