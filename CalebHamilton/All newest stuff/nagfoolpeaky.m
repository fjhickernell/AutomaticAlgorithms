function [out naginfo]=nagfoolpeaky(x)
%This function is being made due to the fact that NAG can't use function
%handles so the peakyfunction has been saved to a file and it can be
%reloaded here so that it has it's own .m file.
load('nagPeaks.mat')
[out,~,~, naginfo]= peakyfunction(x,naginfo);
