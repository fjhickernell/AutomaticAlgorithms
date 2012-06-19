function [y, yprime, ydubprime,info]=peakyfunctionfred(x0,info)
%This function takes a vector x0 of x values and loads specific x values from
%a file. This is used most easily in the script for fooling function.
load(info.filename)
info.sortedX=sort(xsample);
info.sortedX=[info.sortedX((diff(info.sortedX))~=0); max(info.sortedX)];%This makes it chebfun compatible
n=length(info.sortedX)-1;
x=x0(:);
m=length(x);
xx=repmat(x,1,n);
sortedXX=repmat(info.sortedX',m,1);
% yy=1/(size(sortedXX,2)-1)*functionName(xx)+((1-((2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1)))./ ...
%     (sortedXX(:,2:n+1)-sortedXX(:,1:n))).^2).^info.p) ...
%     .*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1)));
yy=((4*(xx-sortedXX(:,1:n)).*(sortedXX(:,2:n+1)-xx)./ ...
    ((sortedXX(:,2:n+1)-sortedXX(:,1:n)).^2)).^info.p) ...
    .*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1)));

%%
%------------------First Derivative---------------------
yyprime=(-4*info.p*(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))).*(1-(2*xx-(sortedXX(:,1:n)...
    +sortedXX(:,2:n+1))).^2./((-sortedXX(:,1:n)+sortedXX(:,2:n+1))).^2).^(-1+info.p)...
    ./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2).*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1)));
%primeNorm=norm(yyprime,inf);

%%
%------------------Second Derivative------------------
yydubprime=(16*info.p*(-1+info.p)*info.p*(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))).^2.*...
    (1-(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))).^2./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2)...
    .^(-2+info.p)./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^4-8*info.p*(1-(2*xx-(sortedXX(:,1:n)...
    +sortedXX(:,2:n+1))).^2./(-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2).^(-1+info.p)./...
    (-sortedXX(:,1:n)+sortedXX(:,2:n+1)).^2).*((xx>=sortedXX(:,1:n))&(xx<=sortedXX(:,2:n+1)));
%dubNorm=norm(yydubprime,inf);

%%
%------------------Output--------------------
%ratio=dubNorm/primeNorm;
yprime=sum(yyprime,2);
yprime=info.sign*reshape(yprime,size(x0));
ydubprime=sum(yydubprime,2);
ydubprime=info.sign*reshape(ydubprime,size(x0));
y=info.RegFunc(x)+info.sign*sum(yy,2);
y=reshape(y,size(x0));
