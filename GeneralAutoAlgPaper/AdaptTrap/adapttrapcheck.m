function [Q,param]=adapttrapcheck(f,param)
%adaptive trapezoidal rule multistage algorithm

param.success=false; %failed to finish within budget
param.nmin=max(param.nmin,ceil((param.tau+1)/2)); %minimum number of trapezoids
ntrap=param.nmin; %present number of trapezoids
xpts=linspace(0,1,ntrap+1); %initial sampling points
fpts=f(xpts); %initial function values
bigdif=fpts(ntrap+1)-fpts(1); %difference quotient across the whole interval
sumf=(fpts(1)+fpts(ntrap+1))/2+sum(fpts(2:ntrap)); %trapezoidal rule sum
iterate=true; %keep iterating
while iterate
    tauok=true; %tau large enough for now
    diff1=diff(fpts); %first difference of function values
    diff2=diff(diff1); %second difference of function values
    absfirstdiff=sum(abs(diff1-bigdif/ntrap)); %estimate 1-norm of first derivative
    abssecdiff=ntrap*sum(abs(diff2)); %estimate total variation of first derivative
    mmult0=2; %next ntrap must be at least this large
    denom=absfirstdiff+abssecdiff/(2*ntrap); %quantity used to for necessary test
    if abssecdiff > denom*param.tau %necessary test for integrand in cone
        param.tau=2*abssecdiff/denom;
        param.nmin=max(param.nmin,ceil((param.tau+1)/2)); %minimum number of trapezoids
        if param.nmin>ntrap; %must increase ntrap
            mmult0=max(ceil(param.nmin/ntrap)); %by at least this factor
            tauok=false; %ntrap is too small compared to tau
        end
    end
    if tauok %ntrap still large enough compared to tau
        %estimate trapezoidal rule error
        errest=param.tau*absfirstdiff/(4*ntrap*(2*ntrap-param.tau));
        if errest <= param.tol %error estimate is small enough
            Q=sumf/ntrap; %answer
            param.success=true; %did not exceed cost budget
            iterate=false; %stop iteration
        elseif ntrap > param.nmax/2 %cannot increase ntrap further
            Q=sumf/ntrap; %answer
            iterate=false; %stop iteration
        end
    end
    if iterate %increase ntrap to hopefully satisfy error criterion 
        ntrapold=ntrap;
        mmult1=ceil(sqrt(param.tau*absfirstdiff/(8*param.tol))/ntrapold);
        mmult2=floor(param.nmax/ntrapold); %cannot increase more than this
        mmult=min(max(mmult0,mmult1),mmult2); %multple of ntrap
        ntrap=ntrapold*mmult; %new number of trapezoids
        xnew=repmat((0:ntrapold-1)/ntrapold,mmult-1,1)...
            +repmat((1:mmult-1)'/ntrap,1,ntrapold); % new nodes
        xnew=xnew(:); %in a long vector
        fnew=f(xnew); %new function values
        sumf=sumf+sum(fnew); %update sum of function values
        fnew=reshape(fnew,mmult-1,ntrapold); %merge new values into old
        fpts=[reshape([fpts(1:ntrapold);fnew],1,ntrap) fpts(ntrapold+1)];
    end
end
param.Q=Q;
param.ntrap=ntrap;
param.errest=errest;
        