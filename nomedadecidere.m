clear all
close all
clc

xF=imread("input\21_peppers256.png");

xF=im2double(im2gray(xF));

[n,m]=size(xF);
n=min(n,m);
xF=xF(1:n,1:n)

m=10;% support PSF 
[PSFtilde,center]=psfGauss([m,m],3);

%PSFtilde=fspecial('Gaussian',m,1);
%H_FT=psf2otf(PSFtilde,[n,n]);

H_FT=psf2otf(PSFtilde,[n,n]);

b = real(ifft2(H_FT.*fft2(xF)));

randn('seed',17)
sigma=0.01;
noise = sigma*randn(n,n);
bb=b+noise;

Dh_FT=psf2otf([1,-1],size(bb));
Dv_FT = psf2otf([1;-1],size(bb));
%rng('default')
x0=bb;

maxit=5000;
alpha_0=1;
epsi=10^(-10);
c=0.5;
tau=0.5;
tol=10^(-7);
z0=x0;

mu_min=0.05;

mu_max=3000;

mu=logspace(log10(mu_min),log10(mu_max),1);

mu=exp(6);

x00=x0;
z00=z0;
b=0;
err=zeros(length(mu),1);
W=err;
g=err;
peak_snr=err;
s_signal=err;
for i=1:length(mu)
    tstart=tic;
    [xx,ff]=nesterovdescentgradient(maxit,x0,mu(i),bb,epsi,H_FT,10^(-20));
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
    b=b+toc(tstart);
end

figure
plot(ff,'Linewidth',1)

xlabel("$i$",'Interpreter','latex')
ylabel("$\mathcal{J}(x_i;\beta)$",'interpreter','latex')

%%


[m,j]=min(err);
[m_G,k]=min(g);
[m_W,l]=min(W);

log(mu(j))
log(mu(k))
log(mu(l))

alpha=1;
beta_0=2;
maxitpcg=100000;
tolpcg=10^(-1);
tol=10^(-3);
tol_r=10^(-5);
maxit1=500;
maxit2=60;

toln=10^(-7);
[betaa,fk1,x1,time1,PSNR1,SSIM1]=gnSTVimplicit(2,maxit1,maxit2,x00,alpha,beta_0,bb,epsi,H_FT,tol,xF,...
maxitpcg,tolpcg,sigma,tol_r,toln,2);

%save("Test confronto\pepperssepsi=10^-3sigma=0.01.mat")

%%

mu_mse=exp(betaa);
[betaa,fk2,x2,time2,PSNR2,SSIM2]=gnSTVimplicit(2,maxit1,maxit2,x00,alpha,beta_0,bb,epsi,H_FT,10^(-4),xF,...
maxitpcg,tolpcg,sigma,tol_r,toln,2);

mu_g=exp(betaa);
[betaa,fk3,x3,time3,PSNR3,SSIM3]=gnSTVimplicit(3,maxit1,maxit2,x00,alpha,beta_0,bb,epsi,H_FT,tol,xF,...
maxitpcg,tolpcg,sigma,tol_r,toln,2);

mu_W=exp(betaa);

%save('Grid for one parameter\sigma=0.01peppers.mat')

%%
clear all
close all
clc

load("Grid for one parameter\sigma=0.05man.mat")

%{
r=x1-xF;
mseopt=1/2*norm(r,'fro')^2;


HxFF=real(ifft2(H_FT.*fft2(x2)));
ResF=HxFF-bb;
gopt=((norm(ResF(:),2).^2-sigma^2*(n^2))^2)/2;

HxFF=real(ifft2(H_FT.*fft2(x3)));
ResF=HxFF-bb;
Resnorm=norm(ResF(:));
ccorrelation2=real(ifft2(fft2(ResF).*conj(fft2(ResF))));
f=ccorrelation2(:)/Resnorm^2;
Wopt=norm(f)^2;
%}
figure

subplot(3,1,1)

loglog(mu,err,'Linewidth',1)
hold on
loglog(mu_mse,1/2*fk1(end),'o','Linewidth',1)
hold on
loglog(mu(j),m,'*','Linewidth',1)
xlabel("e^\beta")
ylabel("MSE(\beta)")
legend('MSE(\beta)','MSE(\beta*) (Gauss-Newton)','MSE(\beta*) (grid-search)')
title('Function MSE(\beta)')


subplot(3,1,2)

loglog(mu,g,'Linewidth',1)
hold on
loglog(mu_g,1/2*fk2(end),'o','Linewidth',1)
hold on
loglog(mu(k),m_G,'*','Linewidth',1)
title('Function g(\beta)')


xlabel("e^\beta")
ylabel("g(\beta)")
legend('g(\beta)','g(\beta*) (Gauss-Newton)','g(\beta*) (grid-search)')


subplot(3,1,3)
loglog(mu,W,'Linewidth',1)
hold on
loglog(mu_W,fk3(end),'o','Linewidth',1)
hold on
loglog(mu(l),m_W,'*','Linewidth',1)
xlabel("e^\beta")
ylabel("W(\beta)")
legend('W(\beta)','W(\beta*) (Gauss-Newton)','W(\beta*) (grid-search)')
title('Function W(\beta)')

figure
subplot(3,1,1)
plot(1/2*fk1,'Linewidth',1)
xlabel("i")
ylabel("MSE(\beta_i)")
title("MSE loss decays along  iteration")

subplot(3,1,2)
plot(1/2*fk2,'Linewidth',1)
xlabel("i")
ylabel("g(\beta_i)")
title("Gaussianity loss decays along iteration")

subplot(3,1,3)
plot(fk3,'Linewidth',1)
xlabel("i")
ylabel("W(\beta_i)")
title("Whiteness loss decays along iteration")

figure
subplot(1,2,1)
imshow2(xF)
title("Original image")

subplot(1,2,2)
imshow2(bb)
title("Observed image")

figure
subplot(1,3,1)
imshow2(x1)
title("Optimal restoration by bilevel MSE")

subplot(1,3,2)
imshow2(x2)
title("Optimal restoration by bilevel Gaussianity")

subplot(1,3,3)
imshow2(x3)
title("Optimal restoration by bilevel Whiteness")

figure

subplot(2,1,1)
loglog(mu,peak_snr,'Linewidth',1)
hold on
loglog(mu_mse,psnr(x1,xF),'o','Linewidth',1)
hold on
loglog(mu_g,psnr(x2,xF),'o','Linewidth',1)
hold on
loglog(mu_W,psnr(x3,xF),'o','Linewidth',1)
title("Function PSNR(\beta)")
xlabel('\beta')
ylabel('PSNR(\beta)')
legend('PSNR','PSNR (MSE)','PSNR (Gaussianity)','PSNR (Whiteness)')
ylim([0,max(peak_snr)+2])

subplot(2,1,2)

loglog(mu,s_signal,'Linewidth',1)
hold on
loglog(mu_mse,ssim(x1,xF),'o','Linewidth',1)
hold on
loglog(mu_g,ssim(x2,xF),'o','Linewidth',1)
hold on
loglog(mu_W,ssim(x3,xF),'o','Linewidth',1)
title("Function SSIM(\beta)")
xlabel('\beta')
ylabel('SSIM(\beta)')
legend('SSIM','SSIM(MSE)','SSIM (Gaussianity)','SSIM (Whiteness)')
ylim([0,max(s_signal)+0.1])


figure
subplot(3,1,1)
plot(time1,'Linewidth',1)
xlabel("i")
ylabel("s")
title("Gauss-Newton Iteration-CPU time for MSE")
subplot(3,1,2)
plot(time2,'Linewidth',1)
xlabel("i")
ylabel("s")
title("Gauss-Newton Iteration-CPU time for g")
subplot(3,1,3)
plot(time3,'Linewidth',1)
xlabel("i")
ylabel("s")
title("Gauss-Newton Iteration-CPU time for W")

figure

subplot(3,1,1)
plot(time1,fk1,'Linewidth',1)
xlabel("s")
ylabel("MSE")
title("Gauss-Newton MSE-CPU time for MSE")
subplot(3,1,2)
plot(time2,fk2,'Linewidth',1)
xlabel("s")
ylabel("g")
title("Gauss-Newton g-CPU time for g")
subplot(3,1,3)
plot(time3,fk3,'Linewidth',1)
xlabel("s")
ylabel("W")
title("Gauss-Newton W-CPU time for W")

figure

subplot(3,1,1)
plot(time1,PSNR1,'Linewidth',1)
xlabel("s")
ylabel("PSNR")
title("Gauss-Newton PSNR-CPU time for MSE")
subplot(3,1,2)
plot(time2,PSNR2,'Linewidth',1)
xlabel("s")
ylabel("PSNR")
title("Gauss-Newton PSNR-CPU time for g")
subplot(3,1,3)
plot(time3,PSNR3,'Linewidth',1)
xlabel("s")
ylabel("PSNR")
title("Gauss-Newton PSNR-CPU time for W")

figure

subplot(3,1,1)
plot(time1,SSIM1,'Linewidth',1)
xlabel("s")
ylabel("SSIM")
title("Gauss-Newton SSIM-CPU time for MSE")
subplot(3,1,2)
plot(time2,SSIM2,'Linewidth',1)
xlabel("s")
ylabel("SSIM")
title("Gauss-Newton SSIM-CPU time for g")
subplot(3,1,3)
plot(time3,SSIM3,'Linewidth',1)
xlabel("s")
ylabel("SSIM")
title("Gauss-Newton SSIM-CPU time for W")