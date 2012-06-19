function y=beginfool(x,info)
%This records all the values of x it is asked for in a file named what you
%like (.txt, .mat, or other suffix must be specially typed)

load(info.filename,'xsample')
xsample=[xsample; x(:)];
save(info.filename,'xsample')
if(info.functionName(x)==0)
    y=zeros(size(x));
else
y=info.functionName(x); %This is the function that you are trying to fit peaks into
end


