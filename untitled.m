clear all
close all
clc

load("datapepperssamenoise.mat")

[q,qq]=find(err==min(err(:)));

maxitpcg=100000;

maxit1=500;

beta_0=[0,0];

[betaa,fk1,x1,time1,PSNR1,SSIM1,ff]=gnSTVimplicit2d(3,maxit1,70,x0,1,beta_0,bb,epsi,...
H_FT,10^(-3),xF,maxitpcg,10^(-2),10^(-6),toln,sigma1,sigma2,sigma_mat,3,M,k,l);


