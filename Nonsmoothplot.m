clear all
close all
clc

xF=imread("star.jpeg");
xF=im2double(im2gray(xF));
[n,m]=size(xF);
n=min(n,m);
xF=xF(1:n,1:n);

m=5;% support PSF           
[PSFtilde,~]=psfGauss([m,m],3);

H_FT=psf2otf(PSFtilde,[n,n]);

b = real(ifft2(H_FT.*fft2(xF)));

randn('seed',17)
sigma=0.05;
noise = sigma*randn(n,n);
bb=b+noise;

Dh_FT=psf2otf([1,-1],size(bb));
Dv_FT = psf2otf([1;-1],size(bb));
x0=zeros(n,n);

maxit=200;
alpha_0=1;
epsi=10^(-3);
c=1;
tau=0.5;
tol=10^(-3);
z0=x0;
beta=0.7;
kk=0.1;
mu_min=5;

mu_max=100;

mu=logspace(log10(mu_min),log10(mu_max),100);

alphaa=0.001;

for i=1:length(mu)
    [xx,ff,xmatrix,alphaa,beta,iter,p]=momentumdescentgradient(maxit,x0,z0,alphaa,beta,mu(i),bb,epsi,Dh_FT,Dv_FT,H_FT,tol);
    r=xx-xF;
    err(i)=1/2*norm(r,'fro')^2;
    
    HxFF=real(ifft2(H_FT.*fft2(xx)));
    ResF=HxFF-bb;
    g(i)=((norm(ResF(:),2).^2-sigma^2*(n^2))^2)/2;
    
    Resnorm=norm(ResF(:));
    ccorrelation2=real(ifft2(fft2(ResF).*conj(fft2(ResF))));
    f=ccorrelation2(:)/Resnorm^2;
    W(i)=norm(f)^2;

    peak_snr(i)=psnr(xx,xF);
    s_signal(i)=ssim(xx,xF);  
    
    x0=xx;
    z0=x0;
end

[m,j]=min(err);
[m_G,k]=min(g);
[m_W,l]=min(W);

log(mu(j))
log(mu(k))
log(mu(l))

c1=0.2;
tau1=0.2;
tol=10^(-6);
z0=x0;

alpha_0=1;
beta_0=3.90;
maxit1=500;
maxit2=200;
tol1=10^(-10);
alpha=1;

[beta_gnerr,fk]=gnSTV(1,Dh_FT,Dv_FT,x0,z0,H_FT,bb,xF,beta_0,maxit1,maxit2,tol,alpha,epsi,alphaa,beta,sigma,n)

mu_err=exp(beta_gnerr);

beta_0=3;

[beta_gngauss,fk2]=gnSTV(2,Dh_FT,Dv_FT,x0,z0,H_FT,bb,xF,beta_0,maxit1,maxit2,tol,alpha,epsi,alphaa,beta,sigma,n)

mu_gauss=exp(beta_gngauss);

[beta_gnwh,fk3]=gnSTV(3,Dh_FT,Dv_FT,x0,z0,H_FT,bb,xF,beta_0,maxit1,maxit2,tol,alpha,epsi,alphaa,beta,sigma,n)

mu_whiteness=exp(beta_gnwh);

[x1]=momentumdescentgradient(maxit,x0,z0,alphaa,beta,mu_err,bb,epsi,Dh_FT,Dv_FT,H_FT,tol);

[x2]=momentumdescentgradient(maxit,x0,z0,alphaa,beta,mu_gauss,bb,epsi,Dh_FT,Dv_FT,H_FT,tol);

[x3]=momentumdescentgradient(maxit,x0,z0,alphaa,beta,mu_whiteness,bb,epsi,Dh_FT,Dv_FT,H_FT,tol);

errr=1/2*fk(end);
gaussf=1/2*fk2(end);
whitef=fk3(end);

figure
loglog(mu,err,'Linewidth',1)
hold on
loglog(mu,g,'Linewidth',1)
hold on
loglog(mu,W,'Linewidth',1)

hold on

loglog(mu(j),m,'*','Linewidth',1)

hold on

loglog(mu(k),m_G,'*','Linewidth',1)

hold on

loglog(mu(l),m_W,'*','Linewidth',1)

hold on

loglog(mu_err,errr,'o','Linewidth',1);

hold on

loglog(mu_gauss,gaussf,'o','Linewidth',1);

hold on

loglog(mu_whiteness,whitef,'o','Linewidth',1);


title('MSE(\mu), g(\mu), W(\mu) for Huber-STV')
legend('MSE(\mu)','g(\mu)','W(\mu)','MSE(\mu*) (grid)','g(\mu*) (grid)',...
    'W(\mu*) (grid)','MSE(\mu*) (unrolling STV Gauss-Newton)','g(\mu*) (unrolling STV Gauss-Newton)',...
    'W(\mu*) (unrolling STV Gauss-Newton)')

xlabel('\mu')

figure
plot(1:length(ff),ff,'-o','Linewidth',1)

title('Momentum gradient descent method+Huber-STV iterations and \mu fixed');


figure
subplot(1,2,1)
imshow(xF)
title('True image')
subplot(1,2,2)
imshow(bb)
title('Observed image')

figure
subplot(1,3,1)
imshow(x1)
title('MGD+Huber-STV image+MSE')

subplot(1,3,2)
imshow(x2)
title('MGD+Huber-STV+gaussianity')

subplot(1,3,3)
imshow(x3)
title('MGD+Huber-STV+RWP')

figure

loglog(mu,peak_snr,'Linewidth',1)
hold on
loglog(mu_err,psnr(x1,xF),'o','Linewidth',1)
hold on
loglog(mu_gauss,psnr(x2,xF),'o','Linewidth',1)
hold on
loglog(mu_whiteness,psnr(x3,xF),'o','Linewidth',1)

title('PSNR(\mu) for Huber-STV')
legend('PSNR','PSNR (MSE)','PSNR (gaussianity)','PSNR (whiteness)','location','southwest')
xlabel('\mu')
ylabel('PSNR')

figure

loglog(mu,s_signal,'Linewidth',1)
hold on
loglog(mu_err,ssim(x1,xF),'o','Linewidth',1)
hold on
loglog(mu_gauss,ssim(x2,xF),'o','Linewidth',1)
hold on
loglog(mu_whiteness,ssim(x3,xF),'o','Linewidth',1)

title('SSIM(\mu) for Huber-STV')
legend('SSIM','SSIM (MSE)','SSIM (gaussianity)','SSIM (whiteness)','location','southwest')
xlabel('\mu')
ylabel('PSNR')

figure
subplot(3,1,1)
plot(1:length(fk),(1/2)*fk,'-o','Linewidth',1)
title('Iteration for GN on STV')
xlabel('i')
ylabel('MSE(\mu_i)')

subplot(3,1,2)
plot(1:length(fk2),(1/2)*fk2,'-o','Linewidth',1)
title('Iteration for GN on STV')
xlabel('i')
ylabel('g(\mu_i)')

subplot(3,1,3)
plot(1:length(fk3),fk3,'-o','Linewidth',1)
title('Iteration for GN on STV')
xlabel('i')
ylabel('W(\mu_i)')
