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

testfun=@(x) 4*sin(x)./exp(x/2);

avec=[1/2 1 2]; na=numel(avec);
bvec=[0.1 0.5 0.9]; nb=numel(bvec);
nvec=[10 20 50]; nn=numel(nvec);
xplot=(-3:0.01:3)';
rmserr=zeros(nn,na,nb);
maxabserr=zeros(nn,na,nb);
condK=zeros(nn,na,nb);
cmat=zeros(nn,na,nb);
dmat=zeros(nn,na,nb);
for iin=1:nn
   n=nvec(iin);
   xdata=(-3:(6/(n-1)):3)';
   ydata=testfun(xdata);
   figure(iin)
   for iia=1:na
      a=avec(iia);
      for iib=1:nb
         b=bvec(iib);
         [K,c,d]=gausshermite(xdata,xdata,a,b);
         cmat(iin,iia,iib)=c;
         dmat(iin,iia,iib)=d;
         condK(iin,iia,iib)=cond(K);
         coef=K\ydata;
         Kplot=gausshermite(xplot,xdata,a,b);
         yplotappx=Kplot*coef;
         yplotexact=testfun(xplot);
         err=abs(yplotexact-yplotappx);
         rmserr(iin,iia,iib)=sqrt(mean(err.*err));
         maxabserr(iin,iia,iib)=max(err);
         subplot(na,nb,(iia-1)*nb+iib)
         plot(xdata,ydata,'r.',xplot,yplotexact,'k-',xplot,yplotappx,'b--')
         %title(['$(a,b,c)=' num2str(a) '$, $b=' num2str(b) '$, cond$(K)=' num2str(condK,3) '$']) 
         text(-3,8,[num2str(a) ','  num2str(b) ','  num2str(c,2) ',' int2str(n)])
%          text(-1,-10,[num2str(condK(iin,iia,iib),3) ',' ...
%             num2str(rmserr(iin,iia,iib),3) ',' ...
%             num2str(maxabserr(iin,iia,iib),3) ])
         axis([-3 3 -10 10])
      end
   end
end
cmat
dmat
condK
rmserr
maxabserr