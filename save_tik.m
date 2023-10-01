clear all
close all
clc

load("datamansamenoise.mat")
maxitpcg=100000;
tolpcg=10^(-3);
tol_r=10^(-5);
toln=10^(-7);
rel_res_sf=3;
tol=10^(-3);

beta0_tik=[0,0];


[betaa_tikerr,fktik_err,xxtik_err,time_tikerr,PSNRtik_err,SSIMtik_err,betaavecerr_tik,ff_tikerr]=tikgnSTVimplicit2d(1,500,60,x0,1,beta0_tik,bb,...
H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l);

betaa_tikerr=[betaa_tikerr(k(1)),betaa_tikerr(l(1))];

[betaa_tikg,fktik_g,xxtik_g,time_tikg,PSNRtik_g,SSIMtik_g,betaavectik_g,ff_tikg]=tikgnSTVimplicit2d(2,500,60,x0,1,beta0_tik,bb,...
H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l);

betaa_tikg=[betaa_tikg(k(1)),betaa_tikg(l(1))];

[betaa_tikW,fktik_W,xxtik_W,time_tikW,PSNRtik_W,SSIMtik_W,betaavectik_W,ff_tikW]=tikgnSTVimplicit2d(3,500,60,x0,1,beta0_tik,bb,...
H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l);

betaa_tikW=[betaa_tikW(k(1)),betaa_tikW(l(1))];

save("tikmanssamenoise.mat")

%%
clear all
close all 
clc

load("datapepperssamenoise.mat")
load("tikpepperssamenoise")

maxitpcg=100000;
tolpcg=10^(-3);
tol_r=10^(-6);
toln=10^(-7);
rel_res_sf=3;
tol=10^(-3);
%k=find(M==1);
%l=find(M==0);
beta01=[0,0];

[betaa_err,fk1,x1,time1,PSNR1,SSIM1,betaavec_err]=gnSTVimplicit2d(1,500,60,x0,1,beta01,bb,...
   epsi,H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l);
[betaa_err1,fk12,x12,time12,PSNR12,SSIM12,betaavec_err2]=gnSTVimplicit2d(1,500,60,x0,1,betaa_tikerr,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l);

[betaa_,fk2,x2,time2,PSNR2,SSIM2,betaavec_g]=gnSTVimplicit2d(2,500,60,x0,1,beta01,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l);
[betaa_2,fk22,x22,time22,PSNR22,SSIM22,betaavec_g2]=gnSTVimplicit2d(2,500,60,x0,1,betaa_tikg,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l);

[betaa,fk3,x3,time3,PSNR3,SSIM3,betaavec]=gnSTVimplicit2d(3,500,60,x0,1,beta01,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l);

[betaa2,fk32,x32,time32,PSNR32,SSIM32,betaavec2]=gnSTVimplicit2d(3,500,60,x0,1,betaa_tikW,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l);

save("bilivelpepperssamenoise.mat")