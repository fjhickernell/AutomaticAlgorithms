function [foolfun,info]=constructquadLyness(info)
%This function takes a vector x0 of x values and loads specific x values from
%a file. This is used most easily in the script for fooling function.
try
    load(info.filename) %load nodes recorded
catch me
    error('No filename given')
end
info.sortedX=sort([info.lower; xsample; info.upper]); %add endpoints
% sort and remove duplicates
info.sortedX=[info.sortedX((diff(info.sortedX))~=0); max(info.sortedX)];%This makes it chebfun compatible

leftnodes=info.sortedX(1:5);
hleft=leftnodes(2)-leftnodes(1);
midnodes=info.sortedX(5:9);
hmid=midnodes(2)-midnodes(1);
rightnodes=info.sortedX(9:13);
hright=rightnodes(2)-rightnodes(1);

%%
a=info.a;
piece=@(x,c) 1./(1 + a^2*(x-c).^2);
integpiece=@(c)(atan(a*(1 - c)) + atan(a*c))/a;
%piece=@(x,c) 1./(exp(a*(x - c)) + exp(-a*(x - c)));
%integpiece=@(c) (-atan(exp(a*(-1 + c))) + atan(exp(a*c)))/a;
ncent=info.ncent;
centers=(0:ncent-1)/(ncent-1);

nconst=5;
A=zeros(nconst,ncent); %matrix containing constraints
b=[1; info.qval; zeros(nconst-2,1)];

for j=1:ncent
    %constraint to make integral to be one
    A(1,j)=integpiece(centers(j));
    %constraint to make the quadature info.qval
    A(2,j)=(hleft*(piece(leftnodes(1),centers(j))...
        +4*piece(leftnodes(3),centers(j))...
        +piece(leftnodes(5),centers(j)))...
        +hmid*(piece(midnodes(1),centers(j))...
        +4*piece(midnodes(3),centers(j))...
        +piece(midnodes(5),centers(j)))...
        +hright*(piece(rightnodes(1),centers(j))...
        +4*piece(rightnodes(3),centers(j))...
        +piece(rightnodes(5),centers(j))))*(2/3);
    %constraint to make the left quadratures match
    A(3,j)=piece(leftnodes(1),centers(j))...
        -4*piece(leftnodes(2),centers(j))...
        +6*piece(leftnodes(3),centers(j))...
        -4*piece(leftnodes(4),centers(j))...
        +piece(leftnodes(5),centers(j));
    %constraint to make the middle quadratures match
    A(4,j)=piece(midnodes(1),centers(j))...
        -4*piece(midnodes(2),centers(j))...
        +6*piece(midnodes(3),centers(j))...
        -4*piece(midnodes(4),centers(j))...
        +piece(midnodes(5),centers(j));
    %constraint to make the middle quadratures match
    A(5,j)=piece(rightnodes(1),centers(j))...
        -4*piece(rightnodes(2),centers(j))...
        +6*piece(rightnodes(3),centers(j))...
        -4*piece(rightnodes(4),centers(j))...
        +piece(rightnodes(5),centers(j));
end
condA=cond(A)
optcoef=A\b;
foolfun=@(x) foolfunmaker(x,optcoef,centers,piece);
end

function y=foolfunmaker(x,coef,centers,piecefun)
ncoef=length(coef);
y=zeros(size(x));
for i=1:ncoef
    y=y+coef(i)*piecefun(x,centers(i));
end
end


