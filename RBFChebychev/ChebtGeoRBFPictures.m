%% Generalized Gauss-Hermite kernel

%% Garbage collection and initialization
format compact %remove blank lines from output
format long e %lots of digits
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',30,'defaulttextfontsize',30) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %latex axis labels
tic

testfun=@(x) 4*sin(3*x)./exp(3*x/2);
a=0.5;
bvec=[0.1 0.5 0.9]; nb=numel(bvec);
nvec=[11 21 51]; nn=numel(nvec);
xplot=(-1:0.01:1)';
rmserr=zeros(nn,nb);
maxabserr=zeros(nn,nb);
condK=zeros(nn,nb);
for iin=1:nn
   n=nvec(iin);
%   xdata=(-1:(2/(n-1)):1)'; %equally spaced
   xdata=cos(pi*(1:-(1/(n-1)):0)'); %Chebychev nodes
   ydata=testfun(xdata);
   for iib=1:nb
      b=bvec(iib);
      K=chebygeo(xdata,xdata,a,b);
      condK(iin,iib)=cond(K);
      coef=K\ydata;
      Kplot=chebygeo(xplot,xdata,a,b);
      yplotappx=Kplot*coef;
      yplotexact=testfun(xplot);
      err=abs(yplotexact-yplotappx);
      rmserr(iin,iib)=sqrt(mean(err.*err));
      maxabserr(iin,iib)=max(err);
      figure(1)
      subplot(nn,nb,(iin-1)*nb+iib)
      plot(xdata,ydata,'r.',xplot,yplotexact,'k-',xplot,yplotappx,'b--')
      text(-1,8,[num2str(b) ',' int2str(n)])
      axis([-1 1 -10 10])
      if iin==1
         figure(2)
         subplot(nb,1,iib)
         plot(xplot,Kplot(:,1:2:n),'b-')
         maxK=max(max(Kplot));
         minK=min(min(Kplot));
         text(-2,-1+ceil(minK),[num2str(b) ',' int2str(n)])
         axis([-1 1 floor(-1+minK) ceil(maxK)])
         %keyboard
      end
   end
end
condK
rmserr
maxabserr