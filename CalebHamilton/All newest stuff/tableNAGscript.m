function naginfo=tableNAGscript(naginfo)
%This script is an altered version of the foolNAGscript to be used
%specially with the makeAtable and with twoScriptGraph. This will also end
%up pruding plots but that is so that it will function for twoScriptGraph.
format long
format compact
naginfo.filename='scriptValues.mat'; %This file name will be passed throughout the script
xsample=[];
save(naginfo.filename,'xsample')

%%
peaks=@(x) peakyfunction(x,naginfo);
save('nagPeaks.mat','peaks','naginfo')
fname='snoopNAG';
switch naginfo.quadtype
    case 'd01ah'
        callautoalg=@(name,lower,upper) d01ah(lower,upper,1e-6,name,int64(100));
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
%%
%-----------------------------Integrals---------------------------
if(naginfo.morepeaks)
    upperbnd=naginfo.sortedX(2:end);
    lowerbnd=naginfo.sortedX(1:end-1);
    withbumps=(upperbnd-lowerbnd)/2*(sqrt(pi)*gamma(naginfo.p+1)/gamma(1.5+naginfo.p));
    failintegral=naginfo.secondattempt;
    naginfo.realintegral=naginfo.c*(sum(withbumps))+naginfo.coefficient/(naginfo.degree+1)*...
        (naginfo.upper^(naginfo.degree+1)-naginfo.lower^(naginfo.degree+1));

else
    largeststart=find(diff(naginfo.sortedX)==max(diff(naginfo.sortedX)));
    lowerbnd=naginfo.sortedX(largeststart);
    upperbnd=naginfo.sortedX(largeststart+1);
    withbumps=(upperbnd-lowerbnd)/2*(sqrt(pi)*gamma(naginfo.p+1)/gamma(1.5+naginfo.p));
    failintegral=naginfo.secondattempt;
    naginfo.realintegral=naginfo.c*(sum(withbumps))+naginfo.coefficient/(naginfo.degree+1)*...
        (naginfo.upper^(naginfo.degree+1)-naginfo.lower^(naginfo.degree+1));
end
%relative error
naginfo.error=abs((naginfo.realintegral-failintegral)/naginfo.realintegral);
%%
%--------------------------------Derivatives------------------------
%===========================Infinity Norm=============================
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
naginfo.primemax=max(primemaxvec);
naginfo.dubmax=max(dubmaxvec);
naginfo.inormratio=naginfo.dubmax/(naginfo.primemax+naginfo.coefficient);
%===============================Two Norm==========================
if(naginfo.morepeaks)
    prime2norm=max(1./sqrt(upperbnd-lowerbnd)*2*naginfo.p*pi^.25*sqrt(gamma(2*naginfo.p-1)/...
      gamma(2*naginfo.p+.5)));
    dub2norm=max(1./sqrt(upperbnd-lowerbnd).^3*2^(5/2)*naginfo.p*sqrt((3*sqrt(pi)*(naginfo.p-1)^2*...
     gamma(naginfo.p*2-3))/(gamma(2*naginfo.p-.5))));
    naginfo.ratio2norm=dub2norm/prime2norm;
else
    u=naginfo.sortedX(naginfo.largeststart+1);
    l=naginfo.sortedX(naginfo.largeststart);
    prime2norm=max(1/sqrt(u-l)*2*naginfo.p*pi^.25*sqrt(gamma(2*naginfo.p-1)/...
    gamma(2*naginfo.p+.5)));
    dub2norm=max(1/sqrt(u-l).^3*2^(5/2)*naginfo.p*sqrt((3*sqrt(pi)*(naginfo.p-1)^2*...
      gamma(naginfo.p*2-3))/(gamma(2*naginfo.p-.5))));
    naginfo.ratio2norm=dub2norm/prime2norm;
end

%%
%-----------------------Plot
load(naginfo.filename)
plot(xx,yplot,'b-',naginfo.sortedX,naginfo.RegFunc(naginfo.sortedX),'ro'),title(naginfo.quadtype),...
    xlabel('x'),ylabel('f(x)')