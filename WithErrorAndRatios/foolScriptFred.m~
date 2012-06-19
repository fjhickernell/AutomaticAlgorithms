%This script will use the beginfool() function to find points called for by
%the function of your choosing (quad, quadgk, chebfun, etc.)
%It then fits peaks in between the points it collects
%Just change this script, and you will get a test of beginfool AND peakyfunction

%This script also gives the ratio of the max derivative values. 
%The last addition was finding the "real" integral so we could get an error
clear all
format long
format compact
info.filename='scriptValues.mat'; %This file name will be passed throughout the script
xsample=[];
save(info.filename,'xsample')

%%
%---------------Function and Bounds----------
%Pick the function you want to put peaks into.
info.RegFunc=@(x) (x-sqrt(3)/2).^2
info.lower=0;
info.upper=1;
info.p=2;
info.sign=-1

%% Ways to call a function
fname='fminbnd';
switch fname
    case 'quadgk'
        callautoalg = @(fun,lower,upper) quadgk(fun,lower,upper);
    case 'quad'
        callautoalg = @(fun,lower,upper) quad(fun,lower,upper);
    case 'chebint'
        callautoalg = @(fun,lower,upper) sum(chebfun(fun,[lower upper]));
    case 'fminbnd'
        callautoalg = @(fun,lower,upper) fminbnd(fun,lower,upper);
end
 

%%
tic
y=@(x) beginfoolfred(x,info);

%The follwing functions should be the one you want to fool
Original=callautoalg(y,info.lower,info.upper)

%%
%-------------Sorting--------------
% load(info.FILENAME); %why not put these three lines in peakyfunction???
% sortedX=sort(xsample);
% save(info.FILENAME,'sortedX');

peaks=@(x) peakyfunctionfred(x,info);


%%
%----------------Making the peaky function and plotting-------
x=info.lower:.005:info.upper; %The domain of the peaky function
[yplot, primeplot, dubplot, info]=peaks(x);

% plot(x,yplot,'r-',x,info.RegFunc(x),'bo') %Here the peaks are plotted over the original function
subplot(2,3,1), plot(x,yplot,'b-',x,info.RegFunc(x),'ro'),title('Peaky Function Overlaid'),...
    legend('Peaky Function','Original','Location','EastOutside')

%plots of Derivatives
subplot(2,3,3), plot(x,primeplot),title('First Derivative') %<<<NOTE: Discontinuous
subplot(2,3,5), plot(x,dubplot),title('Second Derivative')  %<<<NOTE: Discontinuous
break

%%
%The following function needs to be the one you are tricking
inaccurate=callautoalg(peaks,info.lower,info.upper);

%if 'inaccurate' = 'Original' then peakyfunction successfully broke quad, or
%quadgk, etc.
toc


%%
%-------------------Ratio Calculations-----------------
[~,~,~,sortedX]=peaks(0);
[~, yyprimemax]=peaks(sortedX(1:(end-1))+diff(sortedX)./4);
primemax=yyprimemax;
[yy yyprimemax yydub]=peaks(sortedX(1:(end-1)));
dubmax=yydub;

ratio=max(dubmax./primemax)

%%
%--------------Error Calculations--------------

withbumps=0;  %But you should be able to calculate this exactly
for n=1:length(sortedX)-1
    %Change the function your checking ONE MORE TIME
    withbumps=withbumps+callautoalg(peaks,sortedX(n),sortedX(n+1));
end
failintegral=inaccurate
realintegral=withbumps
error=abs((realintegral-failintegral)/realintegral)


