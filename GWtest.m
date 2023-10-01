clear all
close all
clc
xF=imread("input/43_man256.png");
xF=imresize(xF,[150,150]);
xF=im2double(im2gray(xF));
[n,m]=size(xF);
n=min(n,m);
xF=xF(1:n,1:n);

PSFtilde=psfGauss([3,3],3);

%PSF = padarray(PSFtilde,[n-m,n-m],0,'post');
%eigH = fft2(circshift(PSF,1-center));
%diagonalizziamo la matrice di blur H tramite FFT.
H_FT=psf2otf(PSFtilde,[n,n]);

b = real(ifft2(H_FT.*fft2(xF)));

randn('seed',17)
sigma=0.05;
noise = sigma*randn(n,n);
bb=b+noise;

Dh_FT=psf2otf([1,-1],size(bb));
Dv_FT = psf2otf([1;-1],size(bb));

Dhxx=real(ifft2(Dv_FT.*fft2(xF)));

DhT_FT=conj(Dh_FT);
DvT_FT=conj(Dv_FT);
DTD_FT=DhT_FT .* Dh_FT + DvT_FT .* Dv_FT;
HTH_FT=conj(H_FT).*H_FT;

errormin=realmax;
xFnorm=norm(xF,'fro');

bbhat=fft2(bb);
HTbb_FT=conj(H_FT).*bbhat;
mu_min = 0.005;
mu_max =3000;
mu = logspace(log10(mu_min), log10(mu_max),1000);
time_grid=0;
err=zeros(length(mu),1);
W=err;
fg=err;
peaksnr=err;
structsim=err;
for j=1:length(mu) %Studiamo la variazione di Relative Standard Error (RSE).
    tstart=tic;
    sol_FT = HTbb_FT./(HTH_FT + DTD_FT/mu(j));
    xFF = real(ifft2(sol_FT));
    peaksnr(j)=psnr(xF,xFF);
    structsim(j)=ssim(xF,xFF);

    err(j) = (1/2)*norm(xFF(:)-xF(:))^2;
    
    fg(j)=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j),sigma,n);
    W(j)=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,mu(j));    
    a=toc(tstart);
    time_grid=time_grid+a;
end

[errormin,l]=min(err);
mu_err=mu(l);
sol_FT = HTbb_FT./(HTH_FT + DTD_FT/mu_err);
y_opt =real(ifft2(sol_FT));

peaksnr_opt=psnr(xF,y_opt);
structsim_opt=ssim(xF,y_opt);

[m,i]=min(fg);

mu_gauss=mu(i);

[mm,j]=min(W);

mu_W=mu(j);
 
sol_FT = HTbb_FT./(HTH_FT+DTD_FT/mu_gauss);
x_G = real(ifft2(sol_FT));
ee_G=norm(xF-x_G,'fro')/xFnorm
peaksnr_G=psnr(xF,x_G);
structsim_G=ssim(xF,x_G);

sol_FT = HTbb_FT./(HTH_FT+DTD_FT/mu_W);
x_DW = real(ifft2(sol_FT));
ee_DW=norm(xF-x_DW,'fro')/xFnorm
peaksnr_W=psnr(xF,x_DW);
structsim_W=ssim(xF,x_DW);

tol=10^-6;
maxit=100;
beta_0=2;

[beta,fk1,time1]=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,1,sigma,n,bb)

[betaa,fk2,time2]=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,2,sigma,n,bb)

%g_GN=gaussianity(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(betaa),sigma,n);

g_GN=1/2*fk2(end);

R_FT=exp(beta)*conj(H_FT)./(DTD_FT+exp(beta).*HTH_FT);
Rf=real(ifft2(R_FT.*bbhat));
%er_GN=(1/2)*norm(Rf(:)-xF(:))^2;
er_GN=1/2*fk1(end)

[beta_GNW,fk3,time3]=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,3,sigma,n,bb)

%GNW=GRWP(H_FT,bb,HTbb_FT,HTH_FT,DTD_FT,exp(beta_GNW));

GNW=fk3(end);

%%
clear all
close all
clc

load("Tikhonov test 1p\PEPPERS0.05.mat")


figure
subplot(1,2,1)
imshow2(xF)
title("True image")
subplot(1,2,2)
imshow2(bb)
title("Observed image")

figure,
subplot(1,3,1)
imshow2(y_opt)
title("Optimal restoration by bilevel MSE")
subplot(1,3,2)
imshow2(x_G)
title("Optimal restoration by bilevel Gaussianity")
subplot(1,3,3)
imshow2(x_DW)
title("Optimal restoration by bilevel Whiteness")

figure
loglog(mu,fg,'Linewidth',1)
hold on
loglog(mu_gauss,min(fg),'k*','Linewidth',1)
legend("g(\beta)","g(\beta*)")
xlabel("e^\beta")
title("Function g(\beta)")

hold on

loglog(mu,W,'Linewidth',1)
hold on
loglog(mu_W,W(j),'g*','Linewidth',1)

hold on

plot(mu,err,'Linewidth',1)
hold on
loglog(mu_err,errormin,'ro','Linewidth',1)
%hold on
%loglog(mu_gauss,ee_G,'ko','Linewidth',1)
%hold on
%loglog(mu_W,ee_DW,'go','Linewidth',1)

xlabel("\beta")
legend("g(\beta)","g(\beta*)","W(\beta)","W(\beta*)","MSE(\beta)","MSE(\beta*)")
title("Function W(\beta),g(\beta) and MSE(\beta)")

figure

subplot(2,1,1)
loglog(mu,peaksnr,'Linewidth',1)

hold on

loglog(mu_err,peaksnr_opt,'o','Linewidth',1)
hold on
loglog(mu_gauss,peaksnr_G,'o','Linewidth',1)
hold on
loglog(mu_W,peaksnr_W,'o','Linewidth',1)

legend('PSNR(\beta)','PSNR (MSE)','PSNR (Gaussianity)','PSNR (Whiteness)','Location','southwest')

ylabel('PSNR(\beta)')
xlabel('e^\beta')
ylim([0,max(peaksnr)+2])

title('Function PSNR(\beta)')

subplot(2,1,2)

loglog(mu,structsim,'Linewidth',1)

hold on

loglog(mu_err,structsim_opt,'o','Linewidth',1)
hold on
loglog(mu_gauss,structsim_G,'o','Linewidth',1)
hold on
loglog(mu_W,structsim_W,'o','Linewidth',1)

legend('SSIM(\beta)','SSIM (MSE)','SSIM (Gaussianity)','SSIM (Whiteness)','Location','southwest')

ylabel('SSIM(\beta)')
xlabel('e^\beta')

title('Function SSIM(\beta)')
ylim([0,max(structsim)+0.1])

[m,k]=min(err);


figure
loglog(mu,err,'Linewidth',1)
hold on
loglog(exp(beta),er_GN,'o','Linewidth',1)
hold on
loglog(mu,err,'Linewidth',1)
hold on
loglog(mu_err,errormin,'o','Linewidth',1)
xlabel('e^\beta')
legend('\epsilon(\beta)','\MSE(\beta*)(Gauss-Newton)','MSE(\beta)','MSE(\beta*) (grid-search)',...
    'location','northwest')
title('Function \epsilon(\beta) and MSE(\beta)')

figure
subplot(3,1,1)
loglog(mu,err,'Linewidth',1)
hold on
loglog(exp(beta),er_GN,'o','Linewidth',1)
hold on
loglog(mu(k),m,'*','Linewidth',1)
xlabel('e^\beta')
ylabel('MSE(\beta)')
legend('MSE(\beta)','MSE(\beta*) (Gauss-Newton)','MSE(\beta*) (grid-search)')

title('Function MSE(\beta)')



subplot(3,1,2)


loglog(mu,fg,'Linewidth',1)
hold on
loglog(exp(betaa),fk2(end),'o','Linewidth',1)
hold on
loglog(mu(i),min(fg),'*','Linewidth',1)


legend('g(\beta)','g(\beta*) (Gauss-Newton)','g(\beta*) (grid-search)')
xlabel('e^\beta')
ylabel('g(\beta)')
title("Function g(\beta)")

subplot(3,1,3)


loglog(mu,W,'Linewidth',1)

hold on

loglog(exp(beta_GNW),GNW,'o','Linewidth',1)

hold on

loglog(mu_W,min(W),'*','Linewidth',1)
legend('W(\beta)','W(\beta*) (Gauss-Newton)','W(\beta*) (grid-search)')
title("Function W(\beta)")

xlabel('e^\beta')
ylabel('W(\beta)')


figure
subplot(3,1,1)
plot(1/2*fk1,'Linewidth',1)
xlabel("i")
ylabel("MSE(\beta_i)")
title("MSE loss decays along iterations")

subplot(3,1,2)
plot(1/2*fk2,'Linewidth',1)
xlabel("i")
ylabel("g(\beta_i)")
title("Gaussianity loss decays along iterations")

subplot(3,1,3)
plot(fk3,'Linewidth',1)
xlabel("i")
ylabel("W(\beta_i)")
title("Whiteness loss decays along iterations")

figure
subplot(3,1,1)
plot(time1)
subplot(3,1,2)
plot(time2)
subplot(3,1,3)
plot(time3)