%Test out fwht
close all, clear all
format long
format compact

n=2^13; 
x=(0:n-1)/n; 
y=sin(x); 

tic
ywt=fwht(y,[],'sequency'); 
toc
tic
yfft=fft(y); 
%yfft=[yfft(1:n/2); yfft(n:-1:n/2+1)];
yfft=yfft(:)';
toc

sqrt(sum(ywt(1:n/2).^2)), sqrt(sum(ywt.^2))
sqrt(sum(abs(yfft(1:n/2)).^2)), sqrt(sum(abs(yfft).^2))

%loglog(1:n,abs(ywt),'b-',1:n,abs(yfft),'k-','linewidth',2)
loglog(abs(yfft))
figure(gcf)
break

% %% 2-d example
% n=2^13; 
% x=net(sobolset(2),n); 
% y=sin(x(:,1))+(10-x(:,2)).^2; 
% 
% tic
% ywt=fwht(y,[],'sequency'); 
% toc
% tic
% yfft=fft(y)/n; 
% yfft=[yfft(1:n/2); yfft(n:-1:n/2+1)];
% yfft=yfft(:)';
% toc
% 
% sqrt(sum(ywt(1:n/2).^2)), sqrt(sum(ywt.^2))
% sqrt(sum(abs(yfft(1:n/2)).^2)), sqrt(sum(abs(yfft).^2))
% 
% figure
% loglog(1:n,abs(ywt),'b-',1:n,abs(yfft),'k-','linewidth',2)
% figure(gcf)