function [y, yprime, ydubprime,info]=peakyfunction(x0,info)
%This function takes a vector x0 of x values and loads specific x values from
%a file. This is used most easily in the script for fooling function.
try
    load(info.filename) %load nodes recorded
catch
    error('No filename given')
end
info.sortedX=sort([info.lower; xsample; info.upper]); %add endpoints
% sort and remove duplicates
info.sortedX=[info.sortedX((diff(info.sortedX))~=0); max(info.sortedX)];%This makes it chebfun compatible
n=length(info.sortedX)-1; %number of peaks
x=x0(:); %where function is to be evaluated

%% -------------------Function values
m=length(x);
xx=repmat(x,1,n);
sortedXX=repmat(info.sortedX',m,1);
if(~info.morepeaks)
info.largeststart=find(diff(info.sortedX)==max(diff(info.sortedX)));
info.largeststart=info.largeststart(1);
yy=info.c*((4*(xx-sortedXX(:,1:n)).*(sortedXX(:,2:n+1)-xx)./ ...
    ((sortedXX(:,2:n+1)-sortedXX(:,1:n)).^2)).^info.p) ...
    .*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1))).*...
    ((xx>=info.sortedX(info.largeststart))&(xx<=info.sortedX(info.largeststart+1)));
else
    yy=info.c*((4*(xx-sortedXX(:,1:n)).*(sortedXX(:,2:n+1)-xx)./ ...
    ((sortedXX(:,2:n+1)-sortedXX(:,1:n)).^2)).^info.p) ...
    .*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1)));
end
y=info.RegFunc(x0)+info.sign*reshape(sum(yy,2),size(x0));

%% ------------------First Derivative---------------------
if(~info.morepeaks)
    yyprime=info.c*(-4*info.p*(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))).*(1-(2*xx-(sortedXX(:,1:n)...
    +sortedXX(:,2:n+1))).^2./((-sortedXX(:,1:n)+sortedXX(:,2:n+1))).^2).^(-1+info.p)...
    ./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2).*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1))).*...
    ((xx>=info.sortedX(info.largeststart))&(xx<=info.sortedX(info.largeststart+1)));
else
    yyprime=info.c*(-4*info.p*(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))).*(1-(2*xx-(sortedXX(:,1:n)...
    +sortedXX(:,2:n+1))).^2./((-sortedXX(:,1:n)+sortedXX(:,2:n+1))).^2).^(-1+info.p)...
    ./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2).*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1)));
end
  yprime=sum(yyprime,2);
% add the slowly varying function
yprime=info.RegFuncprime(x0)+info.sign*reshape(yprime,size(x0));

%% ------------------Second Derivative------------------
if(~info.morepeaks)
    yydubprime=info.c*(16*info.p*(-1+info.p)*info.p*(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))).^2.*...
    (1-(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))).^2./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2)...
    .^(-2+info.p)./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^4-8*info.p*(1-(2*xx-(sortedXX(:,1:n)...
    +sortedXX(:,2:n+1))).^2./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2).^(-1+info.p)./...
    (-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2).*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1))).*...
    ((xx>=info.sortedX(info.largeststart))&(xx<=info.sortedX(info.largeststart+1)));
else
    yydubprime=info.c*(16*info.p*(-1+info.p)*info.p*(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))).^2.*...
    (1-(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))).^2./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2)...
    .^(-2+info.p)./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^4-8*info.p*(1-(2*xx-(sortedXX(:,1:n)...
    +sortedXX(:,2:n+1))).^2./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2).^(-1+info.p)./...
    (-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2).*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1)));
end
ydubprime=sum(yydubprime,2);
ydubprime=info.RegFuncdubprime(x0)+info.sign*reshape(ydubprime,size(x0));
