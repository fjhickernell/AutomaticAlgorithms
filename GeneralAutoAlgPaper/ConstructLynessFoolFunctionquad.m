%This script will use the find the sample points for quad
%and then construct a function that fools it like Lyness says
%Needs  recordsamplepoints.m
%       constructquadLyness.m

clear all, close all
format long
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
info.filename='scriptValues.mat'; %This file name will be passed throughout the script
xsample=[];
save(info.filename,'xsample')

%% Parameters
info.lower=0; %left endpoint of function domain
info.upper=1;
tol=1e-20; 

%% Find where quad samples
fooler=@(x) recordsamplepoints(x,info);
originalquad=quad(fooler,info.lower,info.upper,tol)

%% Construct fooling function
info.a=10;
info.qval=0.9;
info.ncent=100;
[foolfun,info]=constructquadLyness(info); %foolfun is it

%% Plotting the fooling function
xplot=info.lower:.0025:info.upper; %The domain of the fooling function
yplot=foolfun(xplot);
ymax=max(abs(yplot));
yplot=yplot/ymax;
figure;
pink=[1,0.75,0.75];
%fill([0 xplot 1 0],[0 yplot 0 0],pink,'edgecolor','red','linewidth',1)
nfine=length(info.sortedX);
for j=1:(nfine-1)/2;
    aa=info.sortedX(2*j-1);
    bb=info.sortedX(2*j);
    cc=info.sortedX(2*j+1);
    xparab=aa+(cc-aa).*(0:0.002:1);
    yparab=(xparab-bb).*(xparab-cc)*foolfun(aa)./(ymax*(aa-bb).*(aa-cc)) ...
        +(xparab-aa).*(xparab-cc)*foolfun(bb)./(ymax*(bb-aa).*(bb-cc)) ...
        +(xparab-aa).*(xparab-bb)*foolfun(cc)./(ymax*(cc-aa).*(cc-bb));
%     hsim1=fill([aa xparab cc aa],[0 yparab 0 0],pink,'edgecolor','red','linewidth',3);
%     hold on
end
coarseSimpx=info.sortedX(1:2:end);
ncoarse=length(coarseSimpx);
drkgreen=[0 0.7 0];
for j=1:(ncoarse-1)/2;
    aa=coarseSimpx(2*j-1);
    bb=coarseSimpx(2*j);
    cc=coarseSimpx(2*j+1);
    xparab=aa+(cc-aa).*(0:0.002:1);
    yparab=(xparab-bb).*(xparab-cc)*foolfun(aa)./(ymax*(aa-bb).*(aa-cc)) ...
        +(xparab-aa).*(xparab-cc)*foolfun(bb)./(ymax*(bb-aa).*(bb-cc)) ...
        +(xparab-aa).*(xparab-bb)*foolfun(cc)./(ymax*(cc-aa).*(cc-bb));
%     hsim2=plot([aa xparab cc aa],[0 yparab 0 0],...
%         'LineStyle','-','linewidth',3,'color',drkgreen);
%     hold on
end
% % Color elaborate plot
% h=plot(xplot,yplot,'b-',...
%     info.sortedX,foolfun(info.sortedX)/ymax,'r.',...
%     'markersize',30,'linewidth',3);
% axis([0 1 -0.2 1])
% set(gca,'Xtick',(0:0.2:1),'Ytick',(-0.2:0.2:1))
% hleg=legend(h(1),{'$f$'},'Location','Northwest');
% set(hleg,'Interpreter','latex')
% title('Lyness Fooling Function')
% print -depsc foolquadexample.eps

% Black and white simple plot
h=plot(xplot,yplot,'k-',...
     info.sortedX,foolfun(info.sortedX)/ymax,'k.',...
     'markersize',30,'linewidth',3);
axis([0 1 -0.2 1])
set(gca,'Xtick',(0:0.2:1),'Ytick',(-0.2:0.2:1))
hleg=legend(h,{'$f$','data'},'Location','Northwest');
set(hleg,'Interpreter','latex')
%title('Lyness Fooling Function')
print -depsc foolbwquadexample.eps

%% Try quad out on the fooling function
finalquad=quad(foolfun,info.lower,info.upper,tol)/ymax
finalinteg=1/ymax