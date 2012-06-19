function integral=newtraps(func,lower,upper,n)
%This is a better trapezoidal approximation, utilizing the colon notation
%for vectors. Choosing 'n' includes the two end points.
x=lower+(upper-lower)/(n-1)*(0:n);
yx=func(x);
integral=.5*(upper-lower)/(n-1)*(yx(1)+sum(2*yx(2:n-1))+yx(n));