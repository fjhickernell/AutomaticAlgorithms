%This script will use the beginfool() function to find points called for by
%the function of your choosing (quad, quadgk, chebfun, etc.)
%It then fits peaks in between the points it collects
%Just change this script, and you will get a test of beginfool AND peakyfunction

%This script also gives the ratio of the max derivative values. 
%The last addition was finding the "real" integral so we could get an error
clear all
format long
FILENAME='scriptValues.mat'; %This file name will be passed throughout the script
xsample=[];
save(FILENAME,'xsample')

%%
%---------------Function and Bounds----------
%Pick the function you want to put peaks into.
RegFunc=@(x) 0
lower=0;
upper=1;

%%
tic
y=@(x) beginfool(x,FILENAME,RegFunc);

%The follwing functions should be the one you want to fool
Original=quadgk(y,lower,upper)

%%
%-------------Sorting--------------
load(FILENAME);
sortedX=sort(xsample);
save(FILENAME,'sortedX');

%%
%----------------Making the peaky function and plotting-------
x=lower:.005:upper; %The domain of the peaky function
[yplot primeplot dubplot]=peakyfunction(x,FILENAME,RegFunc);

% plot(x,yplot,'r-',x,RegFunc(x),'bo') %Here the peaks are plotted over the original function
subplot(2,3,1), plot(x,yplot,'b-',x,RegFunc(x),'ro'),title('Peaky Function Overlaid'),...
    legend('Peaky Function','Original','Location','EastOutside')

%plots of Derivatives
subplot(2,3,3), plot(x,primeplot),title('First Derivative') %<<<NOTE: Discontinuous
subplot(2,3,5), plot(x,dubplot),title('Second Derivative')  %<<<NOTE: Discontinuous

%%
%The following function needs to be the one you are tricking
innaccurate=quadgk(@(x) peakyfunction(x,FILENAME,RegFunc),lower,upper);

%if 'inaccurate' = 'Original' then peakyfunction successfully broke quad, or
%quadgk, etc.
toc



%%
%-------------------Ratio Calculations-----------------
[yy yyprimemax]=peakyfunction(sortedX(1:(end-1))+diff(sortedX)./4,FILENAME,RegFunc);
primemax=yyprimemax;
[yy yyprimemax yydub]=peakyfunction(sortedX(1:(end-1)),FILENAME,RegFunc);
dubmax=yydub;

ratio=max(dubmax./primemax)

%%
%--------------Error Calculations--------------

n=length(sortedX)-1;
withbumps=0;
for n=1:length(sortedX)-1
    %Change the function your checking ONE MORE TIME
    withbumps=withbumps+quadgk(@(x) peakyfunction(x,FILENAME,RegFunc),...
    sortedX(n),sortedX(n+1));
end
failintegral=innaccurate
realintegral=withbumps
error=abs((realintegral-failintegral)/realintegral)


