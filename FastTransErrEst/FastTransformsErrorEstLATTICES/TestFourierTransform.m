%Test out nu fft
close all, clear all
format long
format compact

%% We test that our nufft gives us the same result as the Fourier transform
n=2^15;
latticeseq_b2('init0');
xpts=latticeseq_b2(1,n);
nfftr=nufft(2,xpts);
fftr=fft(sort(xpts))/n;
error=norm(nfftr-fftr);
plot(abs(nfftr-fftr))
figure
loglog(1:n,abs(nfftr),'r-',1:n,abs(fftr),'b-','linewidth',2)
hold on
loglog(1:n,abs(nfftr),'r-','linewidth',2)
loglog(1:n,abs(fftr),'b-','linewidth',2)
