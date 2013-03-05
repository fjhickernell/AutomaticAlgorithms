clear all;
close all;
clc

format long;

n=1000000;



a=exp((3*rand(n,1)+1)*log(10));

mu=1./(sqrt(2)*a)+rand(n,1).*(1-sqrt(2)./a);

fp1norm = -(exp(-(a.*(mu-1)).^2) + exp(-(a.*mu).^2)-2);

fpp1norm = 2*a.*(2*sqrt(2/exp(1))-a.*((1-mu).*exp(-(a.*(mu-1))...
    .^2)- a.*mu.*exp(-(a.*mu).^2)));

% fp1norm = -exp(-(a.*(mu-1)).^2) - exp(-(a.*mu).^2);
% 
% fpp1norm = 2*a.^3.*(-((1-mu).*exp(-(a.*(mu-1)).^2)...
%     - mu.*exp(-(a.*mu).^2)));

ratio = fpp1norm./fp1norm;

which = find(ratio(:,1)< 10);

prob1 = length(which)/n;

which = find(ratio(:,1)< 100);

prob2 = length(which)/n;

which = find(ratio(:,1)< 1000);

prob3 = length(which)/n;

display(prob1);
display(prob2);
display(prob3);