%Construct Spiky and Fluky Functions for 
%   Ordinary Trapezoidal Rule

%% Garbage collection
format compact %remove blank lines from output
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %latex axis labels

%16 Subinterval Example
%% Initial Data
n=16; % number of trapezoids
xnodes=(0:n)'/n;
zoom=0.2;
xnodeszoom=xnodes(xnodes<=zoom);
pink=[1,0.75,0.75];
ltgrn=[0.75,1,0.75];

%% Spiky Function
fhump=@(x,a,b) -1 + 60*((x - a).^2).*((b - x).^2)./(b - a)^4;
%fhump=@(x,a,b) 3-4*abs(2*x - (a+b))./(b - a);
fallhumps=@(x,n) fhump(mod(n*x,1),0,1);
fnodes=fallhumps(xnodes,n);
fnodeszoom=fallhumps(xnodeszoom,n);
trapspiken=sum(fnodes.*[1; 2*ones(n-1,1); 1])/(2*n)
trapspikenov2=sum(fnodes(1:2:n+1).*[1; 2*ones(n/2-1,1); 1])/n

%Plot spiky function
nplot=500; 
xplot=(0:nplot)'/nplot;
xplotzoom=zoom*(0:nplot)'/nplot;
fplotzoom=fallhumps(xplotzoom,n);
plot(xplotzoom,fplotzoom,'k-',...
    xnodeszoom,fnodeszoom,'k.')
xlabel('$x$')
ylabel('Spiky Integrand')
axis([0 zoom -1.5 3])
eval(['print -depsc SpikyInteg' int2str(n) 'TrapFigcolor.eps'])

%% Fluky Function
bern2=@(x) 1/6 - x.*(1-x);
bern4=@(x) -1/30 + (x.*(1-x)).^2;
ffluke=@(x,n) 1-(15.*(n.^2))*(bern2(x)+(n.^2)*bern4(x));
fnodes=ffluke(xnodes,n);
trapfluken=sum(fnodes.*[1; 2*ones(n-1,1); 1])/(2*n)
trapflukenov2=sum(fnodes(1:2:n+1).*[1; 2*ones(n/2-1,1); 1])/n

%Plot fluky function
fplot=ffluke(xplot,n);
figure
plot(xplot,fplot,'k-',...
    xnodes,fnodes,'k.')
xlabel('$x$')
ylabel('Fluky Integrand')
axis([0 1 -4e4 4e4])
eval(['print -depsc FlukyInteg' int2str(n) 'TrapFigcolor.eps'])
break

%4 Subinterval Example
%% Initial Data
n=6; % number of trapezoids
xnodes=(0:n)'/n;
zoom=1;
xnodeszoom=xnodes(xnodes<=zoom);

%% Spiky Function
fhump=@(x,a,b) -1 + 60*((x - a).^2).*((b - x).^2)./(b - a)^4;
fallhumps=@(x,n) fhump(mod(n*x,1),0,1);
fnodes=fallhumps(xnodes,n);
fnodeszoom=fallhumps(xnodeszoom,n);
trapspiken=sum(fnodes.*[1; 2*ones(n-1,1); 1])/(2*n);
trapspikenov2=sum(fnodes(1:2:n+1).*[1; 2*ones(n/2-1,1); 1])/n;

%Plot spiky function
nplot=500; 
xplot=(0:nplot)'/nplot;
xplotzoom=zoom*(0:nplot)'/nplot;
fplotzoom=fallhumps(xplotzoom,n);
figure
plot(xplotzoom,fplotzoom,'k-',...
    xnodeszoom,fnodeszoom,'b.',...
    'linewidth',2,'markersize',30)
xlabel('$x$','Interpreter','Latex')
ylabel('Spiky Integrand','Interpreter','Latex')
axis([0 zoom -1.5 3])
print -depsc 'SpikyIntegFigcolor4.eps'

%% Fluky Function
bern2=@(x) 1/6 - x.*(1-x);
bern4=@(x) -1/30 + (x.*(1-x)).^2;
ffluke=@(x,n) 1-(300.*(n.^2))*(bern2(x)+(n.^2)*bern4(x))./n^4;
fnodes=ffluke(xnodes,n);
trapfluken=sum(fnodes.*[1; 2*ones(n-1,1); 1])/(2*n)
trapflukenov2=sum(fnodes(1:2:n+1).*[1; 2*ones(n/2-1,1); 1])/n

%Plot fluky function
fplot=ffluke(xplot,n);
figure
xfill=[reshape([xnodes(1:n) xnodes(1:n) xnodes(2:n+1) xnodes(2:n+1) xnodes(1:n)]',5*n,1); 0];
yfill=[reshape([zeros(n,1) fnodes(1:n) fnodes(2:n+1) zeros(n,1) zeros(n,1)]',5*n,1); 0];
fill(xfill,yfill,pink,'edgecolor','red','linewidth',3);
    hold on
xfill=[reshape([xnodes(1:2:n-1) xnodes(1:2:n-1) xnodes(3:2:n+1) xnodes(3:2:n+1) xnodes(1:2:n-1)]',5*n/2,1); 0];
yfill=[reshape([zeros(n/2,1) fnodes(1:2:n-1) fnodes(3:2:n+1) zeros(n/2,1) zeros(n/2,1)]',5*n/2,1); 0];
fill(xfill,yfill,ltgrn,'edgecolor','green','linewidth',3);

plot(xplot,fplot,'k-',...
    xnodes,fnodes,'b.',...
    xnodes(1:2:n+1),fnodes(1:2:n+1),'g-',...
    'linewidth',3,'markersize',30)
xlabel('$x$','Interpreter','Latex')
ylabel('Fluky Integrand','Interpreter','Latex')
axis([0 1 -7.5 10])
eval(['print -depsc FlukyIntegFigcolor' int2str(n) '.eps'])
