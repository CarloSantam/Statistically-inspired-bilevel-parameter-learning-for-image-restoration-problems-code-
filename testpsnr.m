clear all
close all
clc

load("bilivelcompletepeppers.mat")

[~,~,xstar_err1p,~,PNSR_err1p,SSIM_err1p]=gnSTVimplicit(1,500,30,x0,1,0,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,sigma1,tol_r,toln,rel_res_sf);


[~,~,xstar_g1p,~,PNSR_g1p,SSIM_g1p]=gnSTVimplicit(2,500,30,x0,1,0,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,sigma1,tol_r,toln,rel_res_sf);


[~,~,xstar_w1p,~,PNSR_w1p,SSIM_w1p]=gnSTVimplicit(3,500,30,x0,1,0,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,sigma1,tol_r,toln,rel_res_sf);

%save("peppersoneparameter")

%%
clear all
close all
clc


load("peppersoneparameter")

figure

subplot(1,2,1)

imshow2(bb)

hold on

contour(M,'r')

title("Observed Image")

h=legend("\Omega_1 contour")




subplot(1,2,2)

imshow2(xstar_err1p)

hold on

m=M(n-30,1:n);

[g,ggg]=find(m==1);

a=zeros(n);

a(n-30,ggg)=1;

c=zeros(n,n);

[g,gg]=find(m==0);

c(n-30,gg)=1;

hold on

contour(c,'b')

hold on

contour(a,'g')

title("Deblured peppers with 1 parameter MSE (PSNR="+num2str(PNSR_err1p(end))+","+" SSIM="+num2str(SSIM_err1p(end))+")")

legend("\Omega_2 section","\Omega_1 section")
figure


p=xstar_err1p(n-30,1:n)

plot(1:n,p,'b','linewidth',1)

hold on

p=xstar_err1p(n-30,ggg)

plot(ggg,p,'g','linewidth',1)

legend("1d profile in \Omega_2","1d profile in \Omega_1")



