function info=tableScript(info)
%This script is an altered version of the normal fool script to be used
%specially with the makeAtable and with twoScriptGraph. This will also end
%up pruding plots but that is so that it will function for twoScriptGraph.
format long
format compact
info.filename='scriptValues.mat'; %This file name will be passed throughout the script
xsample=[];
save(info.filename,'xsample')

%% Ways to call a function

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


%%
%----------------------Second Attempt--------------------
info.secondattempt=callautoalg(peaks,info.lower,info.upper);
%if 'secondattempt' = 'Original' then peakyfunction successfully broke quad, or
%quadgk, etc.
%%
%--------------------Working with the Derivatives-------------
info.secondzeros=[info.sortedX(2:end)'; info.sortedX(1:end-1)'; (.5*((info.sortedX(2:end)...
    -info.sortedX(1:end-1))/sqrt(2*info.p-1)+info.sortedX(1:end-1)+...
    info.sortedX(2:end)))'; (.5*((-info.sortedX(2:end)+info.sortedX(1:end-1))/...
    sqrt(2*info.p-1)+info.sortedX(1:end-1)+info.sortedX(2:end)))'];
[~,yp]=peaks(info.secondzeros);
primemaxvec=max(abs(yp));
info.thirdzeros=[info.sortedX(2:end)'; info.sortedX(1:end-1)'; ((info.sortedX(1:end-1)+...
    info.sortedX(2:end))/2)'; (.5*(sqrt(3)*(info.sortedX(2:end)-info.sortedX(1:end-1))/...
    sqrt(2*info.p-1)+info.sortedX(1:end-1)+info.sortedX(2:end)))'; (.5*(sqrt(3)*(-info.sortedX(2:end)...
    +info.sortedX(1:end-1))/sqrt(2*info.p-1)+info.sortedX(1:end-1)+info.sortedX(2:end)))'];
[~,~,ydub]=peaks(info.thirdzeros);
dubmaxvec=max(abs(ydub));
info.primemax=max(primemaxvec);
info.dubmax=max(dubmaxvec);
info.inormratio=info.dubmax/(info.primemax+info.coefficient);

    upperbnd=info.sortedX(2:end);
    lowerbnd=info.sortedX(1:end-1);
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
%--------------Real Integral----------------------
if(info.morepeaks)
    withbumps=(upperbnd-lowerbnd)/2*(sqrt(pi)*gamma(info.p+1)/gamma(1.5+info.p));
    failintegral=info.secondattempt;
    info.realintegral=info.c*(sum(withbumps))+info.coefficient/(info.degree+1)*...
        (info.upper^(info.degree+1)-info.lower^(info.degree+1));
else
    largeststart=find(diff(info.sortedX)==max(diff(info.sortedX)));
    lowerbnd=info.sortedX(largeststart);
    upperbnd=info.sortedX(largeststart+1);
    withbumps=(upperbnd-lowerbnd)/2*(sqrt(pi)*gamma(info.p+1)/gamma(1.5+info.p));
    failintegral=info.secondattempt;
    info.realintegral=info.c*(sum(withbumps))+info.coefficient/(info.degree+1)*...
        (info.upper^(info.degree+1)-info.lower^(info.degree+1));
end

%relative error
info.error=abs((info.realintegral-failintegral)/info.realintegral);
%%
%---------------------Plot
plot(x,yplot,'b-',info.sortedX,info.RegFunc(info.sortedX),'ro'),title(info.fname),...
    xlabel('x'),ylabel('f(x)')

