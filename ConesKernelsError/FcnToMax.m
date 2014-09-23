function fval = FcnToMax(x,xnode,c,gamma,index)
switch index
    case 1
        fval = 0;
        for j = 1:size(xnode)
            fval = fval+c(j)*(x-xnode(j))*exp(-gamma*(x-xnode(j))^2);
        end
        fval = abs(fval);
end
end