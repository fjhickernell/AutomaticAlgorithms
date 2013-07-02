function [Q,param]=nonadapttrap(f,param)

xpts=linspace(0,1,param.ntrap+1);
fpts=f(xpts);
sumf=(fpts(1)+fpts(param.ntrap+1))/2+sum(fpts(2:param.ntrap));
Q=sumf/param.ntrap;
param.Q=Q;
