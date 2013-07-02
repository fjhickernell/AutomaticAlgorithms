%This script will use the find the sample points for quad
%and then construct a function that fools it like Lyness says

clear all, close all
format long
format compact
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
plot(xplot,yplot,'b-',info.sortedX,foolfun(info.sortedX)/ymax,'r.',...
    'markersize',20,'linewidth',2)
axis([0 1 min(yplot) 1])
%title('Lyness Fooling Function')
print -depsc foolquadexample.eps

%% Try quad out on the fooling function
finalquad=quad(foolfun,info.lower,info.upper,tol)/ymax
finalinteg=1/ymax