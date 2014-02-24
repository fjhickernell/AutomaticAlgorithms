function [K,c,d]=gausshermite(t,x,a,b)
nt=size(t,1);
nx=size(x,1);
tminx=bsxfun(@minus,t,x');
t2=repmat(t.*t,1,nx);
x2=repmat(x.*x,1,nt);
c=a*sqrt(b/(1-b^2));
d=0.5*(1-(a.*a)*(1-b)/(1+b));
K=exp(-(c.*c)*tminx.^2 + d*(t2+x2'));
%keyboard