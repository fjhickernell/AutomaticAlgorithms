function y=peakyIntegralfunction(x,nodes)
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
   a=widths(i);
   z=midpts(i);
   b=heights(i);
   xmz=xwh-z;
   xmzma=xmz-a;
   xmzpa=xmz+a;
   y(wh)=b*(4.*a^2 + xmz.^2 + xmzma.*abs(xmzma) - xmzpa.*abs(xmzpa));
end
y=reshape(y,sizex);
