clear all
format long
filename='scriptValues.mat';
% if exist(filename,'file')
%     delete(filename);
% end
xsample=[];
save(filename,'xsample')
y=@(x) beginfool(x,filename);
t=quad(y,0,1);
load(filename);
sortedX=sort(xsample);
save(filename,'sortedX');
x=0:.005:1;
yplot=peakyfunction(x);
plot(x,yplot)
r=quad(@(x) peakyfunction(x),0,1)

break
x1=sortedX(1);
x2=sortedX(2);
between=(x2-x1)/2+x1;
yy=@(x) ((1-(2*(x-between)/(x2-x1)).^2).^2) ...
    .*((x>=x1)&(x<=x2));
x=0:.005:1;
yplot=yy(x);
plot(x,yplot)
r=quad(yy,0,1)


