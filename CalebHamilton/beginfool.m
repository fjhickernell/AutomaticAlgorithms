function y=beginfool(x,filename,functionName)
%This records all the values of x it is asked for in a file named what you
%like (.txt, .mat, or other suffix must be specially typed)

load(filename,'xsample')
xsample=[xsample; x(:)];
save(filename,'xsample')
if(functionName(x)==0)
    y=zeros(size(x));
else
y=functionName(x); %This is the function that you are trying to fit peaks into
end


