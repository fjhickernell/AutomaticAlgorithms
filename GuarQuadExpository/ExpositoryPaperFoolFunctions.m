%This script will use the nodes for quad
%	to construct a functions that fools it like Lyness says
%   Needs foolfunmaker.m

%% Garbage collection and initialization
format compact %remove blank lines from output
clear all %clear all variables
close all %close all figures
format short %set format
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %latex axis labels

%% quad's nodes
xquad=0.13579; %number used by quad to split interval into three parts
xleft=[0 xquad/2 xquad 3*xquad/2 2*xquad];
xctr=[2*xquad 1/4+xquad 1/2 3/4-xquad 1-2*xquad];
xrght=[1-2*xquad 1-3*xquad/2 1-xquad 1-xquad/2 1];
xall=[xleft xctr(2:5) xrght(2:5)]';
nnode=length(xall);

%% Plot a spiky function
fbump=@(x) 4^3*((x.*(1-x)).^3).*((x>=0)&(x<=1)); %one bump
xplot=(0:0.002:1)'; %points to plot
spikyfun=@(x) -1+(35/8)*foolfunmaker(x,@(x,c) fbump((x-c(1))/c(2)),...
    ones(nnode-1,1),[xall(1:nnode-1) diff(xall)]); 
    %spiky function
integralspiky=1;
MATintegralspiky=integral(spikyfun,0,1,'AbsTol',1e-14);
MATquadspiky=quad(spikyfun,0,1,1e-14);
h=plot(xplot,spikyfun(xplot),'k-',xall,-ones(nnode,1),'k.');
axis([0 1 -2.5 4])
set(gca,'Ytick',-2:2:4)
legend(h,{'$f$','data'},'location','southeast')
print -deps 'ExpositoryPaperSpikyquad.eps'
fprintf('  Integral of spiky function = %7.5f\n',integralspiky)
fprintf('By MATLAB''s integral routine = %7.5f\n',MATintegralspiky)
fprintf('   But MATLAB''s quad routine = %7.5f\n\n',MATquadspiky)

%% Plot a fluky function
integval=1; %value of the integral
qval=0.7; %value of quad
a=10; %width of basis functions
piece=@(x,c) 1./(1 + a^2*(x-c).^2); %basis functions
integpiece=@(c)(atan(a*(1 - c)) + atan(a*c))/a; %and their integrals
ncent=100; %number of centers
centers=(0:ncent-1)'/(ncent-1); %and their placements
nconst=5; %number of constraints
A=zeros(nconst,ncent); %matrix containing constraints
b=[1; -1; zeros(nconst-2,1)];
hleft=xall(2)-xall(1);
hright=hleft;
hmid=xall(6)-xall(5);
for j=1:ncent
    %constraint to make integral to be one
    A(1,j)=integpiece(centers(j));
    %constraint to make the quadature info.qval
    A(2,j)=(hleft*(piece(xleft(1),centers(j))...
        +4*piece(xleft(3),centers(j))...
        +piece(xleft(5),centers(j)))...
        +hmid*(piece(xctr(1),centers(j))...
        +4*piece(xctr(3),centers(j))...
        +piece(xctr(5),centers(j)))...
        +hright*(piece(xrght(1),centers(j))...
        +4*piece(xrght(3),centers(j))...
        +piece(xrght(5),centers(j))))*(2/3);
    %constraint to make the left quadratures match
    A(3,j)=piece(xleft(1),centers(j))...
        -4*piece(xleft(2),centers(j))...
        +6*piece(xleft(3),centers(j))...
        -4*piece(xleft(4),centers(j))...
        +piece(xleft(5),centers(j));
    %constraint to make the middle quadratures match
    A(4,j)=piece(xctr(1),centers(j))...
        -4*piece(xctr(2),centers(j))...
        +6*piece(xctr(3),centers(j))...
        -4*piece(xctr(4),centers(j))...
        +piece(xctr(5),centers(j));
    %constraint to make the middle quadratures match
    A(5,j)=piece(xrght(1),centers(j))...
        -4*piece(xrght(2),centers(j))...
        +6*piece(xrght(3),centers(j))...
        -4*piece(xrght(4),centers(j))...
        +piece(xrght(5),centers(j));
end
optcoef=pinv(A)*b;
flukyfun=@(x) foolfunmaker(x,piece,optcoef,centers);
scale=1;
scaledfluky=@(x) flukyfun(x)/scale;
integralfluky=integpiece(centers)'*optcoef/scale;
MATintegralfluky=integral(scaledfluky,0,1,'AbsTol',1e-14);
MATquadfluky=quad(scaledfluky,0,1,1e-14);
figure
h=plot(xplot,scaledfluky(xplot),'k-',xall,scaledfluky(xall),'k.');
axis([0 1 -30 40])
set(gca,'Ytick',-30:10:40)
legend(h,{'$f$','data'},'location','north')
print -deps 'ExpositoryPaperFlukyquad.eps'
fprintf('  Condition number of matrix = %7.5f\n',cond(A))
fprintf('  Integral of fluky function = %7.5f\n',integralfluky)
fprintf('By MATLAB''s integral routine = %7.5f\n',MATintegralfluky)
fprintf('   But MATLAB''s quad routine = %7.5f\n\n',MATquadfluky)
