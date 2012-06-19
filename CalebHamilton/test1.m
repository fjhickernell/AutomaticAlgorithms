sigma=0.01;
upper=1;
lower=0;
n=10000;
t=rand;
f=@(x) exp(-(x-t).^2/(2*sigma^2))/(sqrt(2*pi)*sigma);
tic
newInt=newtraps(f,lower,upper,n)
newEnd=toc
tic
quadInt=quad(f,lower,upper)
quadEnd=toc
tic
gkInt=quadgk(f,lower,upper)
gkEnd=toc
