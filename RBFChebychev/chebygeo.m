function K=chebygeo(t,x,a,b)
ttimesx=bsxfun(@times,t,x');
omt2=1-t.*t;
omx2=1-x.*x;
omt2timesomx2=bsxfun(@times,omt2,omx2');
cosplus=ttimesx - omt2timesomx2;
cosminus=ttimesx + omt2timesomx2;
K=1-a + (a*(1-b))*((cosplus-b)./(1+b^2-2*b*cosplus) ...
   + (cosminus-b)./(1+b^2-2*b*cosminus));
%keyboard