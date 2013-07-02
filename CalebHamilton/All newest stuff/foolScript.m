%This script will use the beginfool() function to find points called for by
%the function of your choosing (quad, quadgk, chebfun, etc.)
%It then fits peaks in between the points it collects
%Just change this script, and you will get a test of beginfool AND peakyfunction

%This script also gives the ratio of the max derivative values. 
%The last addition was finding the "real" integral so we could get an error
clear all
format long
format compact
info.filename='scriptValues.mat'; %This file name will be passed throughout the script
xsample=[];
save(info.filename,'xsample')

%%
%---------------Function and Bounds----------
%Pick the function you want to put peaks into.
info.morepeaks=false;
oneplot=true;
info.degree=1;
info.coefficient=20;
info.RegFunc=... %slowly varying function
    @(x) info.coefficient*x.^info.degree;
info.RegFuncprime=... %and its derivative
    @(x) info.degree.*info.coefficient*x.^(info.degree-1); %slowly varying
info.RegFuncdubprime=... %and its second derivative
    @(x) info.degree.*(info.degree-1).*info.coefficient*x.^(info.degree-2) ;
info.lower=0; %left endpoint of function domain
info.upper=1;
info.p=2;
info.c=1;
info.sign=1;

%% Ways to call a function
info.fname='quadgk';
switch info.fname
    case 'quadgk'
        callautoalg = @(fun,lower,upper) quadgk(fun,lower,upper);
    case 'quad'
        callautoalg = @(fun,lower,upper) quad(fun,lower,upper);
    case 'chebint'
        callautoalg = @(fun,lower,upper) sum(chebfun(fun,[lower upper]));
    case 'fminbnd'
        callautoalg = @(fun,lower,upper) fminbnd(fun,lower,upper);
end
 

%%
y=@(x) beginfool(x,info);
%The follwing functions should be the one you want to fool
original=callautoalg(y,info.lower,info.upper);
peaks=@(x) peakyfunction(x,info);

%%
%----------------Making the peaky function and plotting-------
x=info.lower:.0025:info.upper; %The domain of the peaky function
[yplot, primeplot, dubplot, info]=peaks(x);
if(~oneplot)
    %Here the peaks are plotted over the original function
    subplot(2,3,1), plot(x,yplot,'b-',info.sortedX,info.RegFunc(info.sortedX),'ro'),title('Peaky Function Overlaid'),...
    %legend('Peaky Function','Original','Location','EastOutside')

    %plots of Derivatives
    figure(1)
    subplot(2,3,3), plot(x,primeplot),title('First Derivative') 
    subplot(2,3,5), plot(x,dubplot),title('Second Derivative')  
else
    figure(1)
    plot(x,yplot,'b-',info.sortedX,info.RegFunc(info.sortedX),'ro'),title('Peaky Function Overlaid')
end

%%
%----------------------Second Attempt--------------------
info.secondattempt=callautoalg(peaks,info.lower,info.upper);
%if 'secondattempt' = 'Original' then peakyfunction successfully broke quad, or
%quadgk, etc.
%%
%--------------------Working with the Derivatives-------------
upperbnd=info.sortedX(2:end);
lowerbnd=info.sortedX(1:end-1);
%==============================Infinity Norm====================
info.secondzeros=[upperbnd'; lowerbnd'; (.5*((upperbnd...
    -lowerbnd)/sqrt(2*info.p-1)+lowerbnd+...
    upperbnd))'; (.5*((-upperbnd+lowerbnd)/...
    sqrt(2*info.p-1)+lowerbnd+upperbnd))'];
[~,yp]=peaks(info.secondzeros);
primemaxvec=max(abs(yp));
info.thirdzeros=[upperbnd'; lowerbnd'; ((lowerbnd+...
   upperbnd)/2)'; (.5*(sqrt(3)*(upperbnd-lowerbnd)/...
    sqrt(2*info.p-1)+lowerbnd+upperbnd))'; (.5*(sqrt(3)*(-upperbnd...
    +lowerbnd)/sqrt(2*info.p-1)+lowerbnd+upperbnd))'];
[~,~,ydub]=peaks(info.thirdzeros);
dubmaxvec=max(abs(ydub));
primemax=max(primemaxvec);
dubmax=max(dubmaxvec);
info.inormratio=dubmax/(primemax+info.coefficient);

%================================Two Norm==========================
if(info.morepeaks)
    prime2norm=max(1./sqrt(upperbnd-lowerbnd)*2*info.p*pi^.25*sqrt(gamma(2*info.p-1)/...
      gamma(2*info.p+.5)));
    dub2norm=max(1./sqrt(upperbnd-lowerbnd).^3*2^(5/2)*info.p*sqrt((3*sqrt(pi)*(info.p-1)^2*...
     gamma(info.p*2-3))/(gamma(2*info.p-.5))));
    info.ratio2norm=dub2norm/prime2norm;
else
    u=info.sortedX(info.largeststart+1);
    l=info.sortedX(info.largeststart);
    prime2norm=max(1/sqrt(u-l)*2*info.p*pi^.25*sqrt(gamma(2*info.p-1)/...
    gamma(2*info.p+.5)));
    dub2norm=max(1/sqrt(u-l).^3*2^(5/2)*info.p*sqrt((3*sqrt(pi)*(info.p-1)^2*...
      gamma(info.p*2-3))/(gamma(2*info.p-.5))));
    info.ratio2norm=dub2norm/prime2norm;
end
%%
%----------------------------------Ratio(E)---------------------
E=.001:.005:1;
ratiofunc=64./((info.sortedX(info.largeststart+1)-info.sortedX(info.largeststart))^2*...
    (9*sqrt(3))*(3/2*sqrt(pi)*gamma(info.p+1)/gamma(info.p+3/2)*(1./E-1)+1));
figure(3)
plot(E,ratiofunc),title('R(E)'),xlabel('Relative Error'),ylabel('Infinity Norm Ratio')

%%
%------------------------Our Original For-Loop Method--------------
% bumpyint=0;  %But you should be able to calculate this exactly
% for n=1:length(info.sortedX)-1
%     %Change the function your checking ONE MORE TIME
%     bumpyint=bumpyint+callautoalg(peaks,info.sortedX(n),info.sortedX(n+1));
% end
%bumpyint
%%
%--------------Real Integral----------------------
if  strcmp(info.fname,'fminbnd')
    failmin=peaks(info.secondattempt)
    xmins=(lowerbnd+upperbnd)/2;
    minPeaks=min(peaks(xmins));
   whereitat=xmins(find(peaks(xmins)==minPeaks))
    realmin = minPeaks
else

withbumps=(upperbnd-lowerbnd)/2*(sqrt(pi)*gamma(info.p+1)/gamma(1.5+info.p));
failintegral=info.secondattempt;
info.realintegral=info.c*(sum(withbumps))+info.coefficient/(info.degree+1)*...
    (info.upper^(info.degree+1)-info.lower^(info.degree+1));
%relative error
info.error=abs((info.realintegral-failintegral)/info.realintegral);
end
%%
%----------------------------Display----------------------------
info
