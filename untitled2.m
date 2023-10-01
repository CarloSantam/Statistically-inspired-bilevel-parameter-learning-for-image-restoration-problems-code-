clear all
close all
clc

xF=imread("input\43_man256.png");

xF=imresize_old(xF,[150,150]);

xF=im2double(im2gray(xF));

M = roipoly(xF);

%%

n=min(size(xF));

xF=xF(1:n,1:n);

m=3;% support PSF           
[PSFtilde,~]=psfGauss([m,m],3);

H_FT=psf2otf(PSFtilde,[n,n]);

Dh_FT=psf2otf([1,-1],[n,n]);

Dv_FT=psf2otf([1;-1],[n,n]);

epsi=10^(-3);

b = real(ifft2(H_FT.*fft2(xF)));

sigma1=0.5;

sigma2=0.00001;

k=find(M==1);

l=find(M==0);

randn('seed',17)

A=randn(n);

noise1=sigma1*A.*M;

noise2=sigma2*A.*(ones(n)-M);

bb=b+noise2+noise1;

x0=bb;

[n,~]=size(bb);
mu = ones(n);


sigma_mat=1/sigma1*M+1/sigma2*(ones(size(M))-M);

xstar=x0;
mu1=100000000;

mu2=1;

mu=mu1*M+mu2*(ones(n)-M); 
maxitnesterov=5000;
toln=10^(-7);

[xstar,f]=nesterovdescentgradientmp(maxitnesterov,xstar,mu,bb,epsi,H_FT,toln);

plot(f)