function y=snooperIntegral(x,info)
%This records all the values of x it is asked for in a file named what you
%like (.txt, .mat, or other suffix must be specially typed)

load(info.filename,'xsample')
xsample=[xsample; x(:)];
save(info.filename,'xsample')
y=zeros(size(x));


