function integral=trapezoidInt(func,lower,upper,n)
%This approximates the integral of a function from a lower limit to an
%upper limit. The user is able to choose how many points to use including
%the limits (e.g. if you choose only two points, it will just be the
%limits)
start=tic;
y(1)=func(lower);
 delta=(upper-lower)/(n-1);
 x(1)=lower;
 b=lower+delta;
 i=2;
 x(i)=b;
while b<=upper
    b=b+delta;
    i=i+1;
    x(i)=b;
end
for i=2:n
    traps(i-1)=(x(i)-x(i-1))*(func(x(i))+func(x(i-1)))/2;
end
integral=sum(traps);
timer=toc(start)

