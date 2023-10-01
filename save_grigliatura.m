clear all
close all
clc
%{
xF=imread("input\43_man256.png");
xF=imresize_old(xF,[150,150]);
xF=im2double(im2gray(xF));
M = roipoly(xF);
n=min(size(xF));
xF=xF(1:n,1:n);
m=3;% support PSF
[PSFtilde,~]=psfGauss([m,m],3);
H_FT=psf2otf(PSFtilde,[n,n]);
Dh_FT=psf2otf([1,-1],[n,n]);
Dv_FT=psf2otf([1;-1],[n,n]);
epsi=10^(-3);
b = real(ifft2(H_FT.*fft2(xF)));
sigma1=0.05;
sigma2=0.01;
k=find(M==1);
l=find(M==0);
randn('seed',17)
A=randn(n);
noise1=sigma1*A.*M;
noise2=sigma2*A.*(ones(n)-M);
bb=b+noise2+noise1;
%mu_min=;
%mu_max=1500;
%}
load("datamancomplete.mat")
sigma1=0.01;
sigma2=0.01;
k=find(M==1);
l=find(M==0);
randn('seed',17)
A=randn(n);
noise1=sigma1*A.*M;
noise2=sigma2*A.*(ones(n)-M);
bb=b+noise2+noise1;
mu1=linspace(0.05,4000,500);
mu2=linspace(0.05,4000,500);
[MU1,MU2]=meshgrid(mu1,mu2);
%{
x0=bb;
[n,~]=size(bb);
mu = ones(n);
sigma1=0.05;
sigma2=0.01;
noise1=sigma1*A.*M;
noise2=sigma2*A.*(ones(n)-M);
bb=b+noise1+noise2;
sigma_mat=1/sigma1*M+1/sigma2*(ones(size(M))-M);
%}
x0=bb;
xstar=x0;
sigma_mat=1/sigma1*M+1/sigma2*(ones(size(M))-M);
W=zeros(length(mu2),length(mu1))
nonnormalizedW=zeros(length(mu2),length(mu1))
err=zeros(length(mu2),length(mu1));
g=zeros(length(mu2),length(mu1));
PSNR=zeros(length(mu2),length(mu1));
SSIM=zeros(length(mu2),length(mu1));
maxitnesterov=500;
toln=10^(-7);
t=0;
for i=1:length(mu1)
    for j=1:length(mu2)
        t=t+1
       mu=mu1(i)*M+mu2(j)*(ones(n)-M);
       [xstar,f]=nesterovdescentgradientmp(500,xstar,mu,bb,epsi,H_FT,toln);
       err(j,i)=1/2*norm(xstar-xF,'fro')^2;
       Ress=sigma_mat.*real((ifft2(H_FT.*fft2(xstar))))-sigma_mat.*bb;
       ccorrelation=real(ifft2(fft2(Ress).*conj(fft2(Ress))));
       Resnorm=norm(Ress(:));
       fun=ccorrelation./(Resnorm^2);
       W(j,i)=norm(fun(:))^2;
       Ress1=real(ifft2(H_FT.*fft2(xstar)))-bb;
       ccorrelation1=real(ifft2(fft2(Ress1).*conj(fft2(Ress1))));
       Resnorm1=norm(Ress1(:));
       fun1=ccorrelation1(:)/Resnorm1^2;
       nonnormalizedW(j,i)=norm(fun1(:))^2;
       RESS=Ress1(:);
       ag=abs(RESS(k)).^2-sigma1^2;
       cg=abs(RESS(l)).^2-sigma2^2;
       r1=sum(ag);
       r2=sum(cg);
       g(j,i)=1/2*(r1^2+r2^2);
       PSNR(j,i)=psnr(xstar,xF);
       SSIM(j,i)=ssim(xstar,xF);
    end
end

save("datamansamenoise")

