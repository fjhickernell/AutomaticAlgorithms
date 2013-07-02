function [foolfun,info]=constructchebLyness(info)
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

chebnodes=info.sortedX([1:6 8 9 11]);
extranodes=info.sortedX([7 10]);

%%
a=info.a;
piece=@(x,c) 1./(1 + a^2*(x-c).^2);
%piece=@(x,c) 1./(exp(a*(x - c)) + exp(-a*(x - c)));
ncent=info.ncent;
centers=(0:ncent-1)/(ncent-1);

nnodes=length(info.sortedX);
nconst=1+nnodes;
A=zeros(nconst,ncent); %matrix containing constraints
b=[1; zeros(nconst-1,1)];
xmiss=0.4;

for j=1:ncent
    %constraint to make Chebyshev series miss
    A(1,j)=piece(xmiss,centers(j));
    %constraint to make Chebyshev series match at the extra points
    A(2:nconst,j)=piece(info.sortedX,centers(j));
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

% for j=1:ncent
%     chebpiece=chebfun(@(x) piece(x,centers(j)),[0 1],'length',9); 
%     cj=chebpoly(chebpiece); %Chebyshev coefficients
%     %constraint to make Chebshev series mismatch at 0.4
%     A(1,j)=piece(xmiss,centers(j))-chebpiece(xmiss);
%     %constraint to make Chebshev series match at the extra points
%     A(2:3,j)=piece(extranodes',centers(j))-chebpiece(extranodes');
%     %constraint to make the certain coefficients zero
%     A(4:nconst,j)=cj(whczero)';
% end
% 
% 

