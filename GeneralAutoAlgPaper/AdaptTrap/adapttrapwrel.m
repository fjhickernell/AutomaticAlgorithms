function [Q,param]=adapttrapwrel(f,param)
%adaptive trapezoidal rule multistage algorithm
%   with a relative error criterion

%% Check parameters
if isfield(param,'theta') 
%defines importance of absolute (theta=0) versus relative (theta=1) error
    param.theta=min(max(param.theta,0),1); %must be between 0 and 1
else
    param.theta=0; %default value, pure absolute error
end
if isfield(param,'abstol') %absolute error tolerance
    param.abstol=max(param.abstol,0); %must be non-negative
else
    param.abstol=1e-6; %default value
end
if isfield(param,'reltol') %relative error tolerance
    param.reltol=min(max(param.reltol,0),1); %must be between 0 and 1
else
    param.reltol=1e-6; %default value
end
if isfield(param,'Mcone') %width of the cone, measure of robustness
    param.Mcone=max(param.Mcone,1); %must be greater than 1
else
    param.Mcone=10; %default value
end
if isfield(param,'nmin') %minimun number of trapezoids to start
    param.nmin=max(ceil(param.nmin),param.Mcone+1); 
else
    param.nmin=param.Mcone+1; %default value
end
if isfield(param,'nmax') %maximum number of trapezoids allowed
    param.nmax=max(ceil(param.nmax),1); 
else
    param.nmax=1e8; %default value
end

%% Begin trapezoidal rule
ntrap=param.nmin;
xpts=(0:ntrap)/ntrap; %initial sampling points
%xpts=linspace(0,1,param.nmin+1); %initial sampling points
fpts=f(xpts); %initial function values
sumf=(fpts(1)+fpts(ntrap+1))/2+sum(fpts(2:ntrap)); %trapezoidal rule sum
while ntrap<param.nmax
    Qtemp=sumf/ntrap;
    absfirstdiff=sum(abs(diff(fpts))); %estimate 1-norm of first derivative
    errest=param.Mcone*absfirstdiff/(8*ntrap*(ntrap-param.Mcone));
    efftol=(1-param.theta)*param.abstol ...
            +param.theta*param.reltol*max(errest,abs(Qtemp));
    if errest<=efftol ... %tolerance is met
            || ntrap > param.nmax/2 %cannot increase ntrap further
        break
    else %increase ntrap to hopefully satisfy error criterion
        ntrapold=ntrap;
        mmult=max(2,...
            floor(1.1*(param.Mcone/2+sqrt((param.Mcone/2)^2+...
            param.Mcone*absfirstdiff/(8*efftol)))/ntrapold));
        ntrap=ntrapold*mmult;
        xnew=repmat(xpts(1:ntrapold),mmult-1,1)...
            +repmat((1:mmult-1)'/ntrap,1,ntrapold);
        xpts=[reshape([xpts(1:ntrapold); xnew],1,ntrap) xpts(ntrapold+1)];
        fnew=f(xnew(:));
        sumf=sumf+sum(fnew);
        fnew=reshape(fnew,mmult-1,ntrapold);
        fpts=[reshape([fpts(1:ntrapold); fnew],1,ntrap) fpts(ntrapold+1)];
    end
end
Q=Qtemp-param.theta*param.reltol.*sign(Qtemp)*min(errest,abs(Qtemp));
param.Q=Q;
param.ntrap=ntrap;
param.errest=errest;
        