clear all
close all
clc

% Blur alto e rumore alto
xF=imread("skyscraper.jpg");
xF=im2double(im2gray(xF));
[n,m]=size(xF);
n=min(n,m);
xF=xF(1:n,1:n);

m=10;% support PSF           
[PSFtilde,~]=psfGauss([m,m],3);
H_FT=psf2otf(PSFtilde,[n,n]);

b = real(ifft2(H_FT.*fft2(xF)));

randn('seed',17)
sigma=0.5;
noise = sigma*randn(n,n);
bb=b+noise;

Dh_FT=psf2otf([1,-1],size(bb));
Dv_FT = psf2otf([1;-1],size(bb));
DhT_FT=conj(Dh_FT);
DvT_FT=conj(Dv_FT);
DTD_FT=DhT_FT .* Dh_FT + DvT_FT .* Dv_FT;
HTH_FT=conj(H_FT).*H_FT;

errormin=realmax;
xFnorm=norm(xF,'fro');

bbhat=fft2(bb);
HTbb_FT=conj(H_FT).*bbhat;
mu_min = 0.005;
mu_max =10000;
mu = logspace(log10(mu_min), log10(mu_max),100);

for j=1:length(mu) %Studiamo la variazione di Relative Standard Error (RSE).
    sol_FT = HTbb_FT./(HTH_FT + DTD_FT/mu(j));
    xFF = real(ifft2(sol_FT));
    W(j)=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j));
    fg(j)=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j),sigma,n);
    epsilon(j)=(1/2)*norm(xFF(:)-xF(:))^2;
    peak_snr(j)=psnr(xF,xFF);
    signal_sim(j)=ssim(xF,xFF);
end

[err1,i]=min(epsilon);
[err2,j]=min(fg);
[err3,k]=min(W);

beta_0=0;
tol=10^(-10);
maxit=100;

beta1=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,1,sigma,n,bb);
beta2=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,2,sigma,n,bb);
beta3=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,3,sigma,n,bb);

v3=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(beta3));
v2=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(beta2),sigma,n);
sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta1));
xFF = real(ifft2(sol_FT));
v1=(1/2)*norm(xFF(:)-xF(:))^2;
psnr_err=psnr(xF,xFF);
ssim_err=ssim(xF,xFF);

sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta2));
xFF = real(ifft2(sol_FT));
psnr_gauss=psnr(xF,xFF);
ssim_gauss=ssim(xF,xFF);


sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta3));
xFF = real(ifft2(sol_FT));
psnr_W=psnr(xF,xFF);
ssim_W=ssim(xF,xFF);


figure

subplot(3,1,1)
loglog(mu,epsilon,'Linewidth',1)
hold on
loglog(exp(beta1),v1,'o','Linewidth',1)
hold on
loglog(mu(i),err1,'*','Linewidth',1)
xlabel('\mu')
ylabel('\epsilon(\mu)')

legend('\epsilon(\mu)','\epsilon(\mu*)(GN)','\epsilon(\mu*)(grid)')

title('\epsilon(\mu) with high blur and high noise')


subplot(3,1,2)
loglog(mu,fg,'Linewidth',1)
hold on
loglog(exp(beta2),v2,'o','Linewidth',1)
hold on
loglog(mu(j),err2,'*','Linewidth',1)
xlabel('\mu')
ylabel('g(\mu)')
legend('g(\mu)','g(\mu*)(GN)','g(\mu*)(grid)')
title('g(\mu) with high blur and high noise')

subplot(3,1,3)
loglog(mu,W,'Linewidth',1)
hold on
loglog(exp(beta3),v3,'o','Linewidth',1)
hold on
loglog(mu(k),err3,'*','Linewidth',1)
xlabel('\mu')
ylabel('W(\mu)')

legend('W(\mu)','W(\mu*)(GN)','W(\mu*)(grid)')
title('W(\mu) with high blur and high noise')

figure
subplot(2,1,1)
loglog(mu,peak_snr,'Linewidth',1)
hold on
loglog(exp(beta1),psnr_err,'o','Linewidth',1)
hold on
loglog(exp(beta2),psnr_gauss,'o','Linewidth',1)
hold on
loglog(exp(beta3),psnr_W,'o','Linewidth',1)

legend('PSNR','PSNR(MSE)','PSNR(gaussianity)','PSNR(whiteness)')

title("PSNR for high blur and high noise")

subplot(2,1,2)
loglog(mu,signal_sim,'Linewidth',1)
hold on
loglog(exp(beta1),ssim_err,'o','Linewidth',1)
hold on
loglog(exp(beta2),ssim_gauss,'o','Linewidth',1)
hold on
loglog(exp(beta3),ssim_W,'o','Linewidth',1)

legend('SSIM','SSIM(MSE)','SSIM(gaussianity)','SSIM(whiteness)')

title("SSIM for high blur and high noise")

%%
% Blur basso e rumore alto
xF=imread("skyscraper.jpg");
xF=im2double(im2gray(xF));
[n,m]=size(xF);
n=min(n,m);
xF=xF(1:n,1:n);

m=5;% support PSF           
[PSFtilde,~]=psfGauss([m,m],4);
H_FT=psf2otf(PSFtilde,[n,n]);

b = real(ifft2(H_FT.*fft2(xF)));

randn('seed',17)
sigma=0.5;
noise = sigma*randn(n,n);
bb=b+noise;

Dh_FT=psf2otf([1,-1],size(bb));
Dv_FT = psf2otf([1;-1],size(bb));
DhT_FT=conj(Dh_FT);
DvT_FT=conj(Dv_FT);
DTD_FT=DhT_FT .* Dh_FT + DvT_FT .* Dv_FT;
HTH_FT=conj(H_FT).*H_FT;

errormin=realmax;
xFnorm=norm(xF,'fro');

bbhat=fft2(bb);
HTbb_FT=conj(H_FT).*bbhat;
mu_min = 0.005;
mu_max =10000;
mu = logspace(log10(mu_min), log10(mu_max),100);

for j=1:length(mu) %Studiamo la variazione di Relative Standard Error (RSE).
    sol_FT = HTbb_FT./(HTH_FT + DTD_FT/mu(j));
    xFF = real(ifft2(sol_FT));
    W(j)=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j));
    fg(j)=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j),sigma,n);
    epsilon(j)=(1/2)*norm(xFF(:)-xF(:))^2;
    peak_snr(j)=psnr(xF,xFF);
    signal_sim(j)=ssim(xF,xFF);
end

[err1,i]=min(epsilon);
[err2,j]=min(fg);
[err3,k]=min(W);

beta_0=0;
tol=10^(-10);
maxit=100;

beta1=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,1,sigma,n,bb);
beta2=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,2,sigma,n,bb);
beta3=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,3,sigma,n,bb);

v3=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(beta3));
v2=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(beta2),sigma,n);
sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta1));
xFF = real(ifft2(sol_FT));
v1=(1/2)*norm(xFF(:)-xF(:))^2;
psnr_err=psnr(xF,xFF);
ssim_err=ssim(xF,xFF);

sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta2));
xFF = real(ifft2(sol_FT));
psnr_gauss=psnr(xF,xFF);
ssim_gauss=ssim(xF,xFF);


sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta3));
xFF = real(ifft2(sol_FT));
psnr_W=psnr(xF,xFF);
ssim_W=ssim(xF,xFF);


figure

subplot(3,1,1)
loglog(mu,epsilon,'Linewidth',1)
hold on
loglog(exp(beta1),v1,'o','Linewidth',1)
hold on
loglog(mu(i),err1,'*','Linewidth',1)
xlabel('\mu')
ylabel('\epsilon(\mu)')

legend('\epsilon(\mu)','\epsilon(\mu*)(GN)','\epsilon(\mu*)(grid)')

title('\epsilon(\mu) with low blur and high noise')


subplot(3,1,2)
loglog(mu,fg,'Linewidth',1)
hold on
loglog(exp(beta2),v2,'o','Linewidth',1)
hold on
loglog(mu(j),err2,'*','Linewidth',1)
xlabel('\mu')
ylabel('g(\mu)')
legend('g(\mu)','g(\mu*)(GN)','g(\mu*)(grid)')
title('g(\mu) with low blur and high noise')

subplot(3,1,3)
loglog(mu,W,'Linewidth',1)
hold on
loglog(exp(beta3),v3,'o','Linewidth',1)
hold on
loglog(mu(k),err3,'*','Linewidth',1)
xlabel('\mu')
ylabel('W(\mu)')

legend('W(\mu)','W(\mu*)(GN)','W(\mu*)(grid)')
title('W(\mu) with low blur and high noise')

figure
subplot(2,1,1)
loglog(mu,peak_snr,'Linewidth',1)
hold on
loglog(exp(beta1),psnr_err,'o','Linewidth',1)
hold on
loglog(exp(beta2),psnr_gauss,'o','Linewidth',1)
hold on
loglog(exp(beta3),psnr_W,'o','Linewidth',1)

legend('PSNR','PSNR(MSE)','PSNR(gaussianity)','PSNR(whiteness)')

title("PSNR for low blur and high noise")

subplot(2,1,2)
loglog(mu,signal_sim,'Linewidth',1)
hold on
loglog(exp(beta1),ssim_err,'o','Linewidth',1)
hold on
loglog(exp(beta2),ssim_gauss,'o','Linewidth',1)
hold on
loglog(exp(beta3),ssim_W,'o','Linewidth',1)

legend('SSIM','SSIM(MSE)','SSIM(gaussianity)','SSIM(whiteness)')

title("SSIM for low blur and high noise")
%%

% Blur alto e rumore basso

xF=imread("skyscraper.jpg");
xF=im2double(im2gray(xF));
[n,m]=size(xF);
n=min(n,m);
xF=xF(1:n,1:n);

m=20;% support PSF           
[PSFtilde,~]=psfGauss([m,m],4);
H_FT=psf2otf(PSFtilde,[n,n]);

b = real(ifft2(H_FT.*fft2(xF)));

randn('seed',17)
sigma=0.005;
noise = sigma*randn(n,n);
bb=b+noise;

Dh_FT=psf2otf([1,-1],size(bb));
Dv_FT = psf2otf([1;-1],size(bb));
DhT_FT=conj(Dh_FT);
DvT_FT=conj(Dv_FT);
DTD_FT=DhT_FT .* Dh_FT + DvT_FT .* Dv_FT;
HTH_FT=conj(H_FT).*H_FT;

errormin=realmax;
xFnorm=norm(xF,'fro');

bbhat=fft2(bb);
HTbb_FT=conj(H_FT).*bbhat;
mu_min = 0.005;
mu_max =10000;
mu = logspace(log10(mu_min), log10(mu_max),100);

for j=1:length(mu) %Studiamo la variazione di Relative Standard Error (RSE).
    sol_FT = HTbb_FT./(HTH_FT + DTD_FT/mu(j));
    xFF = real(ifft2(sol_FT));
    W(j)=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j));
    fg(j)=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j),sigma,n);
    epsilon(j)=(1/2)*norm(xFF(:)-xF(:))^2;
    peak_snr(j)=psnr(xF,xFF);
    signal_sim(j)=ssim(xF,xFF);
end

[err1,i]=min(epsilon);
[err2,j]=min(fg);
[err3,k]=min(W);

beta_0=0;
tol=10^(-10);
maxit=100;

beta1=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,1,sigma,n,bb);
beta2=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,2,sigma,n,bb);
beta3=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,3,sigma,n,bb);

v3=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(beta3));
v2=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(beta2),sigma,n);
sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta1));
xFF = real(ifft2(sol_FT));
v1=(1/2)*norm(xFF(:)-xF(:))^2;
psnr_err=psnr(xF,xFF);
ssim_err=ssim(xF,xFF);

sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta2));
xFF = real(ifft2(sol_FT));
psnr_gauss=psnr(xF,xFF);
ssim_gauss=ssim(xF,xFF);


sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta3));
xFF = real(ifft2(sol_FT));
psnr_W=psnr(xF,xFF);
ssim_W=ssim(xF,xFF);


figure

subplot(3,1,1)
loglog(mu,epsilon,'Linewidth',1)
hold on
loglog(exp(beta1),v1,'o','Linewidth',1)
hold on
loglog(mu(i),err1,'*','Linewidth',1)
xlabel('\mu')
ylabel('\epsilon(\mu)')

legend('\epsilon(\mu)','\epsilon(\mu*)(GN)','\epsilon(\mu*)(grid)')

title('\epsilon(\mu) with high blur and low noise')


subplot(3,1,2)
loglog(mu,fg,'Linewidth',1)
hold on
loglog(exp(beta2),v2,'o','Linewidth',1)
hold on
loglog(mu(j),err2,'*','Linewidth',1)
xlabel('\mu')
ylabel('g(\mu)')
legend('g(\mu)','g(\mu*)(GN)','g(\mu*)(grid)','location','southwest')
title('g(\mu) with high blur and low noise')

subplot(3,1,3)
loglog(mu,W,'Linewidth',1)
hold on
loglog(exp(beta3),v3,'o','Linewidth',1)
hold on
loglog(mu(k),err3,'*','Linewidth',1)
xlabel('\mu')
ylabel('W(\mu)')

legend('W(\mu)','W(\mu*)(GN)','W(\mu*)(grid)')
title('W(\mu) with high blur and low noise')

figure
subplot(2,1,1)
loglog(mu,peak_snr,'Linewidth',1)
hold on
loglog(exp(beta1),psnr_err,'o','Linewidth',1)
hold on
loglog(exp(beta2),psnr_gauss,'o','Linewidth',1)
hold on
loglog(exp(beta3),psnr_W,'o','Linewidth',1)

legend('PSNR','PSNR(MSE)','PSNR(gaussianity)','PSNR(whiteness)','location','southwest')

title("PSNR for high blur and low noise")

subplot(2,1,2)
loglog(mu,signal_sim,'Linewidth',1)
hold on
loglog(exp(beta1),ssim_err,'o','Linewidth',1)
hold on
loglog(exp(beta2),ssim_gauss,'o','Linewidth',1)
hold on
loglog(exp(beta3),ssim_W,'o','Linewidth',1)

legend('SSIM','SSIM(MSE)','SSIM(gaussianity)','SSIM(whiteness)','location','southwest')

title("SSIM for high blur and low noise")

%%

% Blur basso e rumore basso


% Blur alto e rumore alto
xF=imread("skyscraper.jpg");
xF=im2double(im2gray(xF));
[n,m]=size(xF);
n=min(n,m);
xF=xF(1:n,1:n);

m=5;% support PSF           
[PSFtilde,~]=psfGauss([m,m],4);
H_FT=psf2otf(PSFtilde,[n,n]);

b = real(ifft2(H_FT.*fft2(xF)));

randn('seed',17)
sigma=0.005;
noise = sigma*randn(n,n);
bb=b+noise;

Dh_FT=psf2otf([1,-1],size(bb));
Dv_FT = psf2otf([1;-1],size(bb));
DhT_FT=conj(Dh_FT);
DvT_FT=conj(Dv_FT);
DTD_FT=DhT_FT .* Dh_FT + DvT_FT .* Dv_FT;
HTH_FT=conj(H_FT).*H_FT;

errormin=realmax;
xFnorm=norm(xF,'fro');

bbhat=fft2(bb);
HTbb_FT=conj(H_FT).*bbhat;
mu_min = 0.005;
mu_max =10000;
mu = logspace(log10(mu_min), log10(mu_max),100);

for j=1:length(mu) %Studiamo la variazione di Relative Standard Error (RSE).
    sol_FT = HTbb_FT./(HTH_FT + DTD_FT/mu(j));
    xFF = real(ifft2(sol_FT));
    W(j)=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j));
    fg(j)=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j),sigma,n);
    epsilon(j)=(1/2)*norm(xFF(:)-xF(:))^2;
    peak_snr(j)=psnr(xF,xFF);
    signal_sim(j)=ssim(xF,xFF);
end

[err1,i]=min(epsilon);
[err2,j]=min(fg);
[err3,k]=min(W);

beta_0=5;
tol=10^(-10);
maxit=100;

beta1=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,1,sigma,n,bb);
beta2=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,2,sigma,n,bb);
beta3=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,3,sigma,n,bb);

v3=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(beta3));
v2=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(beta2),sigma,n);
sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta1));
xFF = real(ifft2(sol_FT));
v1=(1/2)*norm(xFF(:)-xF(:))^2;
psnr_err=psnr(xF,xFF);
ssim_err=ssim(xF,xFF);

sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta2));
xFF = real(ifft2(sol_FT));
psnr_gauss=psnr(xF,xFF);
ssim_gauss=ssim(xF,xFF);


sol_FT = HTbb_FT./(HTH_FT + DTD_FT/exp(beta3));
xFF = real(ifft2(sol_FT));
psnr_W=psnr(xF,xFF);
ssim_W=ssim(xF,xFF);


figure

subplot(3,1,1)
loglog(mu,epsilon,'Linewidth',1)
hold on
loglog(exp(beta1),v1,'o','Linewidth',1)
hold on
loglog(mu(i),err1,'*','Linewidth',1)
xlabel('\mu')
ylabel('\epsilon(\mu)')

legend('\epsilon(\mu)','\epsilon(\mu*)(GN)','\epsilon(\mu*)(grid)')

title('\epsilon(\mu) with low blur and low noise')


subplot(3,1,2)
loglog(mu,fg,'Linewidth',1)
hold on
loglog(exp(beta2),v2,'o','Linewidth',1)
hold on
loglog(mu(j),err2,'*','Linewidth',1)
xlabel('\mu')
ylabel('g(\mu)')
legend('g(\mu)','g(\mu*)(GN)','g(\mu*)(grid)','location','southwest')
title('g(\mu) with low blur and low noise')

subplot(3,1,3)
loglog(mu,W,'Linewidth',1)
hold on
loglog(exp(beta3),v3,'o','Linewidth',1)
hold on
loglog(mu(k),err3,'*','Linewidth',1)
xlabel('\mu')
ylabel('W(\mu)')

legend('W(\mu)','W(\mu*)(GN)','W(\mu*)(grid)')
title('W(\mu) with low blur and low noise')

figure
subplot(2,1,1)
loglog(mu,peak_snr,'Linewidth',1)
hold on
loglog(exp(beta1),psnr_err,'o','Linewidth',1)
hold on
loglog(exp(beta2),psnr_gauss,'o','Linewidth',1)
hold on
loglog(exp(beta3),psnr_W,'o','Linewidth',1)

legend('PSNR','PSNR(MSE)','PSNR(gaussianity)','PSNR(whiteness)','location','southwest')

title("PSNR for low blur and low noise")

subplot(2,1,2)
loglog(mu,signal_sim,'Linewidth',1)
hold on
loglog(exp(beta1),ssim_err,'o','Linewidth',1)
hold on
loglog(exp(beta2),ssim_gauss,'o','Linewidth',1)
hold on
loglog(exp(beta3),ssim_W,'o','Linewidth',1)

legend('SSIM','SSIM(MSE)','SSIM(gaussianity)','SSIM(whiteness)','location','southwest')

title("SSIM for low blur and low noise")