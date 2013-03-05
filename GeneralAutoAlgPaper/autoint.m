function [Q,param]= autoint(f,param)

% param.nmin=max(param.nmin,ceil(param.Mcone)+2); %minimum number of trapezoids
param.nmin=ceil((param.tau+1)/2)+1;
xpts=linspace(0,1,param.nmin+1);
fpts=f(xpts);
sumf=(fpts(1)+fpts(param.nmin+1))/2+sum(fpts(2:param.nmin));
ntrap=param.nmin;
while ntrap<param.nmax
    Gf=sum(abs(diff(fpts)));
%     Gf=sum(abs(fpts(2:end)-fpts(1:end-1)));
%     Gf=max(ntrap*abs(fpts(2:ntrap)-fpts(1:ntrap-1)));
    errest=param.tau*Gf/(8*ntrap*(ntrap-param.tau/2));
    if errest <= param.tol || ntrap > param.nmax
        Q=sumf/(ntrap);
        break
    else
        inflafac=max(ceil(1/ntrap*sqrt(param.tau*Gf/(8*param.tol))),2);
        xnew = repmat(xpts(1:end-1),inflafac-1,1)+repmat((1:inflafac-1)'...
            /(inflafac*ntrap),1,ntrap);
        ynew = f(xnew);
        xnew = [xpts(1:end-1); xnew];
        xpts = [xnew(:); xpts(end)]';
        ynew = [fpts(1:end-1); ynew];
        fpts = [ynew(:); fpts(end)]'; 
        ntrap=ntrap*inflafac;
%         sumf=(fpts(1)+fpts(param.nmin+1))/2+sum(fpts(2:param.nmin));
        sumf=(fpts(1)+fpts(ntrap+1))/2+sum(fpts(2:ntrap));
%         sumf=sumf+sum(ynew);
    end
end
param.Q=Q;
param.ntrap=ntrap;
param.errest=errest;
        