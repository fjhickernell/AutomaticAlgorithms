function [Q,param]=adapttrap(f,param)
%adaptive trapezoidal rule multistage algorithm

param.success=false;
param.nmin=max(param.nmin,ceil((param.tau+1)/2)); %minimum number of trapezoids
xpts=linspace(0,1,param.nmin+1); %initial sampling points
fpts=f(xpts); %initial function values
sumf=(fpts(1)+fpts(param.nmin+1))/2+sum(fpts(2:param.nmin)); %trapezoidal rule sum
ntrap=param.nmin; %present number of trapezoids
while ntrap<param.nmax
    absfirstdiff=sum(abs(diff(fpts))); %estimate 1-norm of first derivative
    errest=param.tau*absfirstdiff/(4*ntrap*(2*ntrap-param.tau));
    if errest <= param.tol %success
        Q=sumf/ntrap; %answer
        param.success=true; %success
        break
    elseif ntrap > param.nmax/2 %cannot increase ntrap further
        Q=sumf/ntrap;
        break
    else %increase ntrap to hopefully satisfy error criterion
        ntrapold=ntrap;
        mmult=max(2,...
            ceil(sqrt(param.tau*absfirstdiff/(8*param.tol))/ntrapold));
        ntrap=ntrapold*mmult;
        xnew=repmat((0:ntrapold-1)/ntrapold,mmult-1,1)...
            +repmat((1:mmult-1)'/ntrap,1,ntrapold);
        xnew=xnew(:);
        fnew=f(xnew);
        sumf=sumf+sum(fnew);
        fnew=reshape(fnew,mmult-1,ntrapold);
        fpts=[reshape([fpts(1:ntrapold);fnew],1,ntrap) fpts(ntrapold+1)];
    end
end
param.Q=Q;
param.ntrap=ntrap;
param.errest=errest;
        