function TrianglePeak(bwcolor)
%Plot one and two triangle peak functions

%% Garbage collection
format compact %remove blank lines from output
format long %lots of digits
close all %close all figures
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLegendInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %latex axis labels

%% Initial Data
tic
if strcmp(bwcolor,'bw') %For black and white
    linecolor='k';
else %For color
    linecolor='b';
end

%16 Subinterval Example
%% Initial Data
a=0;
b=1;

%% One Peak
t=0.25;
h=0.2;
disp('One Peak')
peak=@(x,t,h) (x>=t & x<=t+2*h).*(h-abs(x-t-h));
nplot=500; 
xplot=a+(0:nplot)'*((b-a)/nplot);
fplot=peak(xplot,t,h);
plot(xplot,fplot,[linecolor '-'])
xlabel('$x$')
ylabel(['$\mbox{peak}(x,' num2str(t) ',' num2str(h) ')$'])
axis([0 1 -0.05 h+0.05])
print('-depsc',['OnePeakFig' bwcolor '.eps'])

%% Two Peaks
hcut=0.2;
t=0.65;
h=0.1;
b=3/4;
disp('TwoPeaks')
twopeak=@(x,t,h,up) peak(x,0,hcut)+up*b*peak(x,t,h);
fplotup=twopeak(xplot,t,h,1);
figure
plot(xplot,fplotup,[linecolor '-'])
xlabel('$x$')
ylabel(['$\mbox{twopk}(x,' num2str(t) ',' num2str(h) ',+)$'])
axis([0 1 -0.05 hcut+0.05])
print('-depsc',['TwoPeakFig' bwcolor '.eps'])
toc



