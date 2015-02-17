function GaussProbTrapExample(bwcolor)
%Integrate Gaussian probability density by the
%   trapezoidal rule

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
disp('Gaussian probability density example')
gray=0.9* [1 1 1];
ltgray=0.5*[1 1 1];
ltblue=[0.4 0.4 1];
ltltblue=[0.85 0.85 1];
if strcmp(bwcolor,'bw') %For black and white
    linecolor='k'; pointcolor='k'; dashcolor=ltgray; shadecolor=gray;
else %For color
    linecolor='b'; pointcolor='r'; dashcolor=ltblue; shadecolor=ltltblue; 
end
n=4; % number of trapezoids
xnodes=(0:n)'/n;

%% Gaussian Probability or erf
sig=1/2;
fGauss=@(x) exp(-x.*x/(2*sig^2))/(sqrt(2*pi)*sig);
fpGauss=@(x) -x.*exp(-x.*x/(2*sig^2))/(sqrt(2*pi)*sig.^3);
fppGauss=@(x) (-1+x.*x/(sig.^2)).*exp(-x.*x/(2*sig^2))/(sqrt(2*pi)*sig.^3);
fnodes=fGauss(xnodes);
trapfGaussn=sum(fnodes.*[1; 2*ones(n-1,1); 1])/(2*n)
exactinteg=2*normcdf(1/sig)-1
err=exactinteg-trapfGaussn
varfprime=-2*fpGauss(sig)+fpGauss(1)
normfpp=integral(@(x) abs(fppGauss(x)),0,1)
errupbd=varfprime/(8*n*n)
tolcoef=sqrt(varfprime/8)
guar_integ=integral_g(fGauss)

%% Plot Gaussian probability function
nplot=500; 
xplot=(0:nplot)'/nplot;
fplot=fGauss(xplot);
figure
xfill=[xnodes; xnodes(n+1:-1:1)];
yfill=[fnodes; zeros(n+1,1)];
h=zeros(n+6,1);
h(1)=fill(xfill,yfill,shadecolor,'edgecolor',shadecolor);
hold on
h(2:3)=plot(xplot,fplot,[linecolor '-'],...
    xnodes,fnodes,[pointcolor '.']);
h(4:n+6)=plot([xnodes xnodes]',[zeros(n+1,1) fnodes]','--', ...
    [0 1],[0 0],'--',xnodes,fnodes,'--','color',dashcolor);
xlabel('$x$')
%ylabel('$f_{\mbox{easy}}(x)$')
set(gca,'xtick',0:0.25:1)
axis([0 1 0 0.9])
legend(h([2 1]),{'$f_{\mbox{easy}}$',['$T_{' int2str(n) '}(f_{\mbox{easy}})$']},...
    'location',[0.6 0.65 0.35 0.25])
legend('boxoff')
eval(['print -depsc GaussInteg' int2str(n) 'TrapFig' bwcolor '.eps'])

%% Heuristic Algorithm
tol=1e-4;
n=2;
fudge=1.2;
xnodes=(0:1/n:1);
fnodes=fGauss(xnodes);
trapnov2=(fnodes(1)+fnodes(3))/2;
trapn=(fnodes(1)+2*fnodes(2)+fnodes(3))/4;
errest=fudge*(trapn-trapnov2)/3;
while errest>tol;
   oldn=n;
   n=2*n;
   xnew=xnodes(1:oldn)+1/n;
   fnew=fGauss(xnew);
   tempnodes=[xnodes(1:oldn);xnew];
   xnodes=[tempnodes(:)' xnodes(oldn+1)];
   tempf=[fnodes(1:oldn);fnew];
   fnodes=[tempf(:)' fnodes(oldn+1)];
   trapnov2=trapn;
   trapn=((fnodes(1)+fnodes(n+1))/2+sum(fnodes(2:n)))/n;
   errest=fudge*(trapn-trapnov2)/3;
end
trapn
errguess=abs(exactinteg-trapn)
toc
   

