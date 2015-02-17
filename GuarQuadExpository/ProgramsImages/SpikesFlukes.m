function SpikesFlukes(bwcolor)
%Construct Spiky and Fluky Functions for 
%   Ordinary Trapezoidal Rule

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
    linecolor='k'; pointcolor='k'; 
else %For color
    linecolor='b'; pointcolor='r'; 
end

%16 Subinterval Example
n=16; % number of trapezoids
xnodes=(0:n)'/n;
zoom=1;
xnodeszoom=xnodes(xnodes<=zoom);

%% Spiky Function
disp('Spiky example for trapezoidal rule')
fhump=@(x,a,b) (15/8)*n*((x - a).^2).*((b - x).^2)./(b - a)^4; %max value is (15/8)*n/16
fallhumps=@(x,n) fhump(mod(n*x,1),0,1);
fnodes=fallhumps(xnodes,n);
fnodeszoom=fallhumps(xnodeszoom,n);
trapspiken=sum(fnodes.*[1; 2*ones(n-1,1); 1])/(2*n)
trapspikenov2=sum(fnodes(1:2:n+1).*[1; 2*ones(n/2-1,1); 1])/n
errestspike=abs(trapspiken-trapspikenov2)/3
integralspike=integral(@(x) fallhumps(x,n),0,1)
%guarspike1=integral_g(@(x) fallhumps(x,n),'ninit',n+1)
%guarspike2=integral_g(@(x) fallhumps(x,n),'ninit',2*n+1)

%Plot spiky function
nplot=500; 
xplot=(0:nplot)'/nplot;
xplotzoom=zoom*(0:nplot)'/nplot;
fplotzoom=fallhumps(xplotzoom,n);
plot(xplotzoom,fplotzoom,[linecolor '-'],...
    xnodeszoom,fnodeszoom,[pointcolor '.'])
xlabel('$x$')
ylabel(['$f_{\mbox{spiky}}(x;' int2str(n) ')$'])
axis([0 zoom -0.1 2])
eval(['print -depsc SpikyInteg' int2str(n) 'TrapFig' bwcolor '.eps'])

%% Fluky Function
disp('Fluky example for trapezoidal rule')
bern2=@(x) 1/6 - x.*(1-x);
bern4=@(x) -1/30 + (x.*(1-x)).^2;
ffluke=@(x,n) 1-((15/2).*(n.^2))*(bern2(x)+(n.^2)*bern4(x));
fnodes=ffluke(xnodes,n);
trapfluken=sum(fnodes.*[1; 2*ones(n-1,1); 1])/(2*n)
trapflukenov2=sum(fnodes(1:2:n+1).*[1; 2*ones(n/2-1,1); 1])/n
errestfluke=abs(trapfluken-trapflukenov2)/3
guarfluke=integral_g(@(x) ffluke(x,n),0,1,1e-7)

%Plot fluky function
fplot=ffluke(xplot,n);
figure
%plot(xplot,fplot,[linecolor '-'],xnodes,fnodes,[pointcolor '.'])
plot(xplot,fplot,[linecolor '-'])
xlabel('$x$')
ylabel(['$f_{\mbox{fluky}}(x;' int2str(n) ')$'])
%ylabel('Fluky Integrand')
axis([0 1 -2e4 2e4])
eval(['print -depsc FlukyInteg' int2str(n) 'TrapFig' bwcolor '.eps'])

%% Big Not Fluky Function
disp('Big example for trapezoidal rule')
fbig=@(x,n) 1 -((15/2).*(n.^4))*bern4(x);
fnodes=fbig(xnodes,n);
trapbign=sum(fnodes.*[1; 2*ones(n-1,1); 1])/(2*n)
trapbignov2=sum(fnodes(1:2:n+1).*[1; 2*ones(n/2-1,1); 1])/n
errestbig=abs(trapbign-trapbignov2)/3
guarbig=integral_g(@(x) fbig(x,n),0,1,1e-7)

%Plot big function
fplot=fbig(xplot,n);
figure
%plot(xplot,fplot,[linecolor '-'],xnodes,fnodes,[pointcolor '.'])
plot(xplot,fplot,[linecolor '-'])
xlabel('$x$')
ylabel(['$f_{\mbox{big}}(x;' int2str(n) ')$'])
axis([0 1 -2e4 2e4])
eval(['print -depsc BigInteg' int2str(n) 'TrapFig' bwcolor '.eps'])
toc

