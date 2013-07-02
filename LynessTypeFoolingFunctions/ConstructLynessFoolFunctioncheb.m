%This script will use the find the sample points for cheb
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

%% Find where cheb samples
fooler=@(x) recordsamplepoints(x,info);
originalquad=chebfun(fooler,[info.lower info.upper]);

%% Construct fooling function
info.a=5;
info.ncent=12;
[foolfun,info]=constructchebLyness(info); %foolfun is it

%% Try cheb out on the fooling function
finalcheb=chebfun(foolfun,[info.lower info.upper])

%% Plotting the fooling and the original function
xplot=info.lower:.0025:info.upper; %The domain of the fooling function
yplot=foolfun(xplot);
figure;
plot(xplot,yplot,'b-',xplot,finalcheb(xplot),'k--',...
    info.sortedX,foolfun(info.sortedX),'r.',...
    'markersize',20,'linewidth',2)
title('Lyness Fooling Function')

