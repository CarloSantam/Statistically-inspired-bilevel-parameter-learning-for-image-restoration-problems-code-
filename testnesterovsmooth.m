clear all
close all
clc

xF=imread("cameraman.jpg");
xF=im2double(im2gray(xF));
[n,m]=size(xF);
n=min(n,m);
xF=xF(1:n,1:n);

m=10 %support PSF           
[PSFtilde,~]=psfGauss([m,m],3);
%PSF = padarray(PSFtilde,[n-m,n-m],0,'post');
%eigH = fft2(circshift(PSF,1-center));
%diagonalizziamo la matrice di blur H tramite FFT.
H_FT=psf2otf(PSFtilde,[n,n]);

b = real(ifft2(H_FT.*fft2(xF)));

randn('seed',17)
sigma=0.03;
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

errormin=realmax;
xFnorm=norm(xF,'fro');

bbhat=fft2(bb);
HTbb_FT=conj(H_FT).*bbhat;
mu=36;
sol_FT = HTbb_FT./(HTH_FT + DTD_FT/mu);
xFF = real(ifft2(sol_FT));
c=0.4;
tau=0.4;
tol=10^(-5);
maxit=70;
rng('default')
x0=rand(n,n);
y0=x0;
alpha_0=1;
[x_n,ff,kk]=smoothnesterovdescentgradient(maxit,x0,y0,alpha_0,mu,bb,Dh_FT,Dv_FT,H_FT,c,tau,tol,xFF);

z=real(ifft2(H_FT.*fft2(xFF)))-bb;
Dhx=real(ifft2(Dh_FT.*fft2(xFF)));
Dvx=real(ifft2(Dv_FT.*fft2(xFF)));

f=mu*norm(z,'fro')^2+norm(Dhx,'fro')^2+norm(Dvx,'fro')^2;
figure
loglog(1:length(kk),kk/norm(xFF,'fro'),'Linewidth',1)
title('NGD error for l2-l2 Tikhonov')

figure
loglog(1:length(ff),ff,'-o','Linewidth',1)
title('Iteration for l2-l2 Tikhonov + Huber descent gradient')


mu_min=0.05;
mu_max=200;

mu=logspace(log10(mu_min),log10(mu_max),1);

for i=1:length(mu)
   xFF=smoothnesterovdescentgradient(maxit,x0,y0,alpha_0,mu(i),bb,Dh_FT,Dv_FT,H_FT,c,tau,tol,xFF);
   err(i)=1/2*norm(xF-xFF,'fro')^2;
   xFF = real(ifft2(HTbb_FT./(HTH_FT + DTD_FT/mu(i))));
   err_2(i)=1/2*norm(xF-xFF,'fro')^2;
end
figure
loglog(mu,err,'Linewidth',1)
hold on
loglog(mu,err_2,'Linewidth',1)

[m,i]=min(err);
log(mu(i))

xlabel('\mu')
ylabel('MSE(\mu)')

legend('Nesterov MSE','"classic" MSE')

hold on

c=0.4;
tau=0.4;
c1=0.1;
tau1=0.1;
tol=10^(-3);
rng('default')
x0=rand(n,n);
tol1=10^(-10);

y0=x0;

alpha_0=1;
beta_0=0;
maxit1=30;
maxit2=70;

[beta,fk]=gnsmoothtest(Dh_FT,Dv_FT,x0,y0,H_FT,bb,xF,xFF,beta_0,maxit1,maxit2,tol,tol1,c,tau,c1,tau1,alpha_0)

xFF=smoothnesterovdescentgradient(maxit,x0,y0,alpha_0,exp(beta),bb,Dh_FT,Dv_FT,H_FT,c,tau,tol,xFF);

errr=1/2*norm(xF-xFF,'fro')^2;

loglog(exp(beta),errr,'o','Linewidth',1)

legend('Nesterov MSE','"classic" MSE','MSE(\mu*)(GN + NGD)')

figure

plot(1:length(fk),fk,'-o','Linewidth',1)
title('Gauss Newton + NGD (MSE) Iteration')
xlabel('\mu')
ylabel('MSE(\mu_i)')

figure
subplot(1,3,1)
imshow(xF)
title('Real image')
subplot(1,3,2)
imshow(bb)
title('Observed image')
subplot(1,3,3)
imshow(x_n)
title('l2-l2 TIK restored image with MSE')