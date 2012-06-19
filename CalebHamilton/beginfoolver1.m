function y=beginfool(x,filename)
%This records all the values of x it is asked for in a file named what you
%like (.txt, .mat, or other suffix must be specially typed)
% fileID=fopen(filename,'a+');
% fprintf(fileID,'%f\n ',x);
load(filename,'xsample')
xsample=[xsample; x(:)];
save(filename,'xsample')
s=size(x);
y=zeros(s);


