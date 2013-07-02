function [fspline,param]=adaptlinspline(f,param)
%adaptive trapezoidal rule multistage algorithm

param.success=false;
param.nmin=max(param.nmin,ceil((param.tau+1)/2)); %minimum number of trapezoids
xpts=linspace(0,1,param.nmin+1); %initial sampling points
fpts=f(xpts); %initial function values
nint=param.nmin; %present number of intervals
while nint<param.nmax
    absfirstdiff=nint*max(abs(diff(fpts))); %estimate inf-norm of first derivative
    errest=param.tau*absfirstdiff/(4*nint*(2*nint-param.tau));
    if errest <= param.tol %success
        fspline=@(x) ppval(interp1((0:nint)/nint,fpts,'linear','pp'),x);
        param.success=true; %success
        break
    elseif nint > param.nmax/2 %cannot increase nint further
        fspline=@(x) ppval(interp1((0:nint)/nint,fpts,'linear','pp'),x);
        break
    else %increase nint to hopefully satisfy error criterion
        nintold=nint;
        mmult=max(2,...
            ceil(sqrt(param.tau*absfirstdiff/(8*param.tol))/nintold));
        nint=nintold*mmult;
        xnew=repmat((0:nintold-1)/nintold,mmult-1,1)...
            +repmat((1:mmult-1)'/nint,1,nintold);
        xnew=xnew(:);
        fnew=f(xnew);
        fnew=reshape(fnew,mmult-1,nintold);
        fpts=[reshape([fpts(1:nintold);fnew],1,nint) fpts(nintold+1)];
    end
end
param.fspline=fspline;
param.nint=nint;
param.errest=errest;
        