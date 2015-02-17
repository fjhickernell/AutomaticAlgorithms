function FoolIntegral(bwcolor)
%Fool MATLAB's integral.m

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
if strcmp(bwcolor,'bw') %For black and white
    linecolor='k'; pointcolor='k'; 
else %For color
    linecolor='b'; pointcolor='r'; 
end

%% Nice Example
tic
integral(@(x) exp(-x.*x),0,1)

%% Set up database file
info.filename='scriptValuesquadgk.mat'; %This file name will be passed throughout the script
xsample=[];
save(info.filename,'xsample')

%% Call the automatic algorithm to fool
Original=integral(@(x) snooperIntegral(x,info),0,1,...
   'AbsTol',1e-14,'RelTol',1e-14);
load(info.filename,'xsample')

%% Gauss-Konrod nodes and weights taken from quadgk
pnodes = [ ...
    0.2077849550078985; 0.4058451513773972; 0.5860872354676911; ...
    0.7415311855993944; 0.8648644233597691; 0.9491079123427585; ...
    0.9914553711208126];
pwt = [ ...
    0.2044329400752989, 0.1903505780647854, 0.1690047266392679, ...
    0.1406532597155259, 0.1047900103222502, 0.06309209262997855, ...
    0.02293532201052922];
pwt7 = [0,0.3818300505051189,0,0.2797053914892767,0,0.1294849661688697,0];
allnodes = [-pnodes(end:-1:1); 0; pnodes];
allwt = [pwt(end:-1:1), 0.2094821410847278, pwt];
errwt = allwt - [pwt7(end:-1:1), 0.4179591836734694, pwt7];

%% Variable transform
hscale=@(t,a,b) 0.25*(b-a)*t.*(3-t.^2)+0.5*(b+a);
hp=@(t,a,b) 0.75*(b-a)*(1-t.^2);
numint=10;

%% Construct a peaky integrand
disp('Spiky integrand for integral.m')
peakf=@(x) spikyIntegralfunction(x,xsample);
figure
xplot=(0:2e-6:0.01)';
fplot=peakf(xplot);
% h=plot(xplot,fplot,'b-',...
%    xplot([1 end]),zeros(2,1),'k-',...
%    xsample,peakf(xsample),'r.');
h=plot(xplot,fplot,[linecolor '-'],...
   xsample,peakf(xsample),[pointcolor '.']);
axis([0 0.01 -0.2 2])
set(gca,'Ytick',(0:0.4:2))
%xlabel('$x$')
%ylabel('$f_{\mbox{peaky}}$')
%hleg=legend(h,{'$f_{\mbox{spiky}}$','data'},'Location','Southeast');
%legend(gca,'boxoff');
%set(hleg,'Interpreter','latex')
print('-depsc',['SpikyFoolIntegral' bwcolor '.eps'])

peakintegral=integral(peakf,0,1,'AbsTol',1e-13,'RelTol',1e-13)

%% Prepare all nodes and weights
initnodes=bsxfun(@plus,0.1*allnodes,(-0.9:0.2:0.9));
initnodes=initnodes(:);
initnodesh=hscale(initnodes,0,1);

hpsample=hp(initnodes,0,1);
initwts=repmat(allwt,1,10);
initwts=initwts'.*hpsample;
initerrwts=repmat(errwt,1,10);
initerrwts=initerrwts'.*hpsample;

%% Construct a fluky integrand
disp('Fluky integrand for integral.m')
a=40; %peakiness of pieces
piece=@(x,c) 1./(1 + a^2*(x-c).^2);
integpiece=@(c)(atan(a*(1 - c)) + atan(a*c))/a;
ncent=50;
centers=(0:ncent-1)/(ncent-1);
nconst=12;
A=zeros(nconst,ncent); %matrix containing constraints
b=[1; 0.9999; zeros(nconst-2,1)];
for j=1:ncent
   %constraint to make integral to be one
   A(1,j)=integpiece(centers(j));
   %constraint to make the quadature something else
   piecenodes=piece(xsample,centers(j));
   A(2,j)=sum(piecenodes.*initwts)/10;
   for jj=1:10
      %constraint to make the error estimate zero
      wh=(jj-1)*15+(1:15);
      A(2+jj,j)=sum(piecenodes(wh).*initerrwts(wh));
   end
end
condA=cond(A)
optcoef=pinv(A)*b;
flukef=@(x) flukyIntegralfunction(x,piece,centers,optcoef);
xplot=(0:0.002:1);
fscale=max(abs(flukef(xplot)));
optcoef=optcoef/fscale;

integf=integpiece(centers)*optcoef

flukef=@(x) flukyIntegralfunction(x,piece,centers,optcoef);
figure
%h=plot(xplot,flukef(xplot),'b-',xsample,flukef(xsample),'r.');
h=plot(xplot,flukef(xplot),[linecolor '-'],...
    xsample,flukef(xsample),[pointcolor '.']);
axis([0 1 -1 1.2])
set(gca,'Ytick',(-1:0.5:1))
%xlabel('$x$')
%ylabel('$f_{\mbox{peaky}}$')
% hleg=legend(h,{'$f_{\mbox{fluky}}$','data'},'Location','Southeast');
% legend(gca,'boxoff');
% set(hleg,'Interpreter','latex')
print('-depsc',['FlukyFoolIntegral' bwcolor '.eps'])

flukeintegral=quadgkFredLook(flukef,0,1,'AbsTol',1e-13,'RelTol',1e-13)

toc

