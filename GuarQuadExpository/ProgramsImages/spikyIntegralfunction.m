function y=spikyIntegralfunction(x,nodes)
%This function takes a vector of nodes 
%  and creates a bumpy function

sizex=size(x);
nnodes=numel(nodes);
midpts=0.5*(nodes(2:nnodes)+nodes(1:nnodes-1));
widths=diff(nodes)/4;
heights=1./(2.*widths.*widths);
nx=numel(x);
y=zeros(nx,1);
for i=1:nnodes-1
   wh=(x>nodes(i))&(x<nodes(i+1));
   xwh=x(wh);
   a=nodes(i);
   b=nodes(i+1);
   y(wh)=(30/((b-a).^4))*((xwh-a).*(b-xwh)).^2;
end
y=reshape(y,sizex);
