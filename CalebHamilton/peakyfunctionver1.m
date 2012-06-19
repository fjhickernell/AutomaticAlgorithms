function y=peakyfunction(x0)
load scriptValues.mat
n=length(sortedX)-1;
x=x0(:);
m=length(x);
xx=repmat(x,1,n);
sortedXX=repmat(sortedX',m,1);
keyboard
yy=((1-(2*xx-(sortedXX(:,1:n)+sortedXX(:,2:n+1))./ ...
    (sortedXX(:,2:n+1)-sortedXX(:,1:n))).^2).^2) ...
    .*((x>=sortedXX(:,1:n))&(x<=sortedXX(:,2:n+1)));
y=sum(yy,2);
y=reshape(y,size(x0));
