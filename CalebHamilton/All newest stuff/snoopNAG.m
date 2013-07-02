function out=snoopNAG(x)
load('nagPeaks.mat')
load(naginfo.filename,'xsample')
xsample=[xsample; x(:)];
save(naginfo.filename,'xsample')
out=naginfo.RegFunc(x);