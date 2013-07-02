%This script runs multiple algorithms multiple times with a change in the
%parameters each time. Choose 'd' or 'p' as your change parameter. Eitehr
%the degree of the input function or the p value of the peak functions will
%be incremented by one. You can also choose how many times you want to run
%the algorithms by changing the 'times' variable. The larger the 'times'
%value, the longer it will take.

clear all
tic
tableinfo.E=.001:.005:1;
tableinfo.morepeaks=false;
tableinfo.lower=0;
tableinfo.upper=1;
tableinfo.c=1;
tableinfo.p=2;
tableinfo.coefficient=20;
tableinfo.degree=1;
tableinfo.sign=1;
tableinfo.RegFunc=@(x) tableinfo.coefficient.*x.^tableinfo.degree;
tableinfo.RegFuncprime=... %and its derivative
    @(x) tableinfo.degree.*tableinfo.coefficient*x.^(tableinfo.degree-1); %slowly varying
tableinfo.RegFuncdubprime=... %and its second derivative
    @(x) tableinfo.degree.*(tableinfo.degree-1).*tableinfo.coefficient*x.^(tableinfo.degree-2);

change='e'; %What parameter will you increment?
times=10;  %How many times will you iterate?
figure(1)
i=1;
ratiovecah=[];
for(counter=1:times)
 
    subplot(2,3,1)
    hold on
    tableinfo.quadtype='d01ah';
    naginfo=tableNAGscript(tableinfo);
    xah(i)=naginfo.error;
    yah(i)=naginfo.inormratio;
    hold off

    subplot(2,3,2)
    hold on
    tableinfo.quadtype='d01aj';
    naginfo=tableNAGscript(tableinfo);
    xaj(i)=naginfo.error;
    yaj(i)=naginfo.inormratio;
    hold off

    subplot(2,3,3)
    hold on
    tableinfo.quadtype='d01ak';
    naginfo=tableNAGscript(tableinfo);
    xak(i)=naginfo.error;
    yak(i)=naginfo.inormratio;
    hold off

    subplot(2,3,4)
    hold on
    tableinfo.fname='quadgk';
    info=tableScript(tableinfo);
    xgk(i)=info.error;
    ygk(i)=info.inormratio;
    hold off

    subplot(2,3,5)
    hold on
    tableinfo.fname='quad';
    info=tableScript(tableinfo);
    xqu(i)=info.error;
    yqu(i)=info.inormratio;
    hold off

    subplot(2,3,6)
    hold on
    tableinfo.fname='chebint';
    info=tableScript(tableinfo);
    xch(i)=info.error;
    ych(i)=info.inormratio;
    hold off
    
    i=i+1;
switch change
    case 'p'
        tableinfo.p=tableinfo.p+1;
    case 'e'
        tableinfo.coefficient=tableinfo.coefficient*2;
    case 'c'
        tableinfo.c=tableinfo.c+5;
    case 'd'
        tableinfo.degree=tableinfo.degree+1;
end
end

%%
%-------------------------------------Plot---------------------------------
% splineah=@(t)spline(xah,yah,t); 
% splineaj=@(t)spline(xaj,yaj,t);
% splineak=@(t)spline(xak,yak,t);
% splinequ=@(t)spline(xqu,yqu,t);
% splinegk=@(t)spline(xgk,ygk,t);
% splinech=@(t)spline(xch,ych,t);
% t=0:.01:1;
figure(2)
hold on
% semilogy(t,splineah(t),'r-',t,splineaj(t),'b-',t,splineak(t),'g-',...
%     t,splinequ(t),'k-',t,splinegk(t),'c-',t,splinech(t),'y-')
semilogy(xah,yah,'ro',xaj,yaj,'bo',xak,yak,'go',xgk,ygk,'co',xqu,yqu,'ko',xch,ych,'yo'),...
    legend('d01ah','d01aj','d01ak','quadgk','quad','chebfun','Location','Best',...
    'linewidth','MarkerFaceColor','red',4)
% ylim([0 1e8]) 
xlabel('Relative error');
ylabel('Inf-Norm');
hold off
toc