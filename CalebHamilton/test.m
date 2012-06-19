%This test was trying to replicate the examples in Prof. Hickernell's
%power point "Monte Carlo Algorithms Where the Integrand is Unknown"

% 'f' represents the final function on slide 4/28. I can't replicate
%his answers using quad()......
clear all
f=@(x) 1+cos(200*pi*x);
g=@(x) (2/sqrt(pi)).*exp(-x.^2);
start1=tic;
q=quad(f,0,1);
end1=toc(start1);
start2=tic;
t=quad(g,0,1);
end2=toc(start2);
end1
end2
format long
q
t
