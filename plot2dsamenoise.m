
close all
clear all
clc


G="bilivelmansamenoise.mat";

load(G)

%load("datapepperscomplete(grigliapiùfine).mat")

[minimo_err]=min(err(:));

[mu2_caerr, mu1_caerr] = find(err == minimo_err);

figure
contourf(log(MU1),log(MU2),log(err),'LineWidth',1,'EdgeColor','none')
hold on
scatter(log(mu1(mu1_caerr)),log(mu2(mu2_caerr)),'r o','LineWidth',1)
hold on
scatter(betaa_err(k(1)),betaa_err(l(1)),'*','LineWidth',1)
hold on
scatter(betaa_err1(k(1)),betaa_err1(l(1)),'+','LineWidth',1)
hold on
plot(betaavec_err(1,:),betaavec_err(2,:),'k .-','LineWidth',1)
hold on
plot(betaavec_err2(1,:),betaavec_err2(2,:),'b.-','LineWidth',1)

legend("MSE($\mu_1$,$\mu_2$)","MSE($\mu_1^*$,$\mu_2^*$) (grid)","MSE($\mu_1^*$,$\mu_2^*$) (GN)","MSE($\mu_1^*$,$\mu_2^*$) (GN TIK +TV)","GN trajectory","GN (TIK+TV) trajectory",'Location','southwest','Interpreter','Latex')

xlabel("log($\mu_1$)",'Interpreter','Latex')
ylabel("log($\mu_2$)",'Interpreter','Latex')

title("MSE level curve")


figure

surfc(log(MU1),log(MU2),log(err),'LineWidth',1,'EdgeColor', 'none')
%hold on
%contour3(log(MU1),log(MU2),log(err),'LineWidth',1)
hold on
plot3(log(mu1(mu1_caerr)),log(mu2(mu2_caerr)),log(err(mu2_caerr,mu1_caerr)),'r o','LineWidth',1)
hold on
plot3(betaa_err(k(1)),betaa_err(l(1)),log(1/2*fk1(end)),'g *','LineWidth',1)
hold on
plot3(betaa_err1(k(1)),betaa_err1(l(1)),log(1/2*fk12(end)),'g ^','LineWidth',1)
hold on
plot3((betaavec_err(1,:)),(betaavec_err(2,:)),log(1/2*fk1),'k .-','LineWidth',1)
hold on
plot3(betaavec_err2(1,:),betaavec_err2(2,:),log(1/2*fk12),'b .-','LineWidth',1)



xlabel("log($\mu_1$)",'Interpreter','Latex')
ylabel("log($\mu_2$)",'Interpreter','Latex')
zlabel("log(MSE($\mu_1$,$\mu_2$))",'Interpreter','Latex')

title("MSE 3d plot")

legend("MSE($\mu_1$,$\mu_2$)","MSE($\mu_1^*$,$\mu_2^*$) level curve","MSE($\mu_1^*$,$\mu_2^*$) (grid)","MSE($\mu_1^*$,$\mu_2^*$) (GN)","MSE($\mu_1^*$,$\mu_2^*$) (GN TIK+TV)","GN trajectory","GN (TIK+TV) trajectory",'Interpreter','Latex')

[minimo_g]=min(g(:));

[mu2_cag, mu1_cag] = find(g == minimo_g);

figure
contourf(log(MU1),log(MU2),log(g),'LineWidth',1,'EdgeColor','none')
hold on
scatter(log(mu1(mu1_cag)),log(mu2(mu2_cag)),'r o','LineWidth',1)
hold on
scatter(betaa_(k(1)),betaa_(l(1)),'*','Linewidth',1)
hold on
scatter(betaa_2(k(1)),betaa_2(l(1)),'g o','Linewidth',1)
hold on
plot(betaavec_g(1,:),betaavec_g(2,:),'k .-','LineWidth',1)
hold on
plot(betaavec_g2(1,:),betaavec_g2(2,:),'b .-','LineWidth',1)


legend("g($\mu_1$,$\mu_2$)","g($\mu_1^*$,$\mu_2^*$) (grid)","g($\mu_1^*$,$\mu_2^*$) (GN)","g($\mu_1^*$,$\mu_2^*$) (GN TIK+TV)","GN trajectory","GN trajectory (GN TIK+TV)",'Interpreter','Latex')

xlabel("log($\mu_1$)",'Interpreter','Latex')
ylabel("log($\mu_2$)",'Interpreter','Latex')

title("Gaussianity level curve")

figure

surfc(log(MU1),log(MU2),log(g),'LineWidth',1,'EdgeColor', 'none')
%hold on
%contour3(log(MU1),log(MU2),log(g),'LineWidth',1)
hold on
plot3(log(mu1(mu1_cag)),log(mu2(mu2_cag)),log(g(mu2_cag,mu1_cag)),'r o','LineWidth',1)
hold on
plot3((betaa_(k(1))),(betaa_(l(1))),log(1/2*fk2(end)),'*','LineWidth',1)
hold on
plot3((betaa_2(k(1))),(betaa_2(l(1))),log(1/2*fk22(end)),'g o','LineWidth',1)
hold on
plot3((betaavec_g(1,:)),(betaavec_g(2,:)),log(1/2*fk2),'k .- ','LineWidth',1)
hold on
plot3((betaavec_g2(1,:)),(betaavec_g2(2,:)),log(1/2*fk22),'b .- ','LineWidth',1)


xlabel("log($\mu_1$)",'Interpreter','Latex')
ylabel("log($\mu_2$)",'Interpreter','Latex')
zlabel("log(g($\mu_1$,$\mu_2$))",'Interpreter','Latex')

legend("g($\mu_1$,$\mu_2$)","g($\mu_1^*$,$\mu_2^*$) level curve","g($\mu_1^*$,$\mu_2^*$) (grid)","g($\mu_1^*$,$\mu_2^*$) (GN)","g($\mu_1^*$,$\mu_2^*$) (GN TIK + TV)","GN trajectory","GN trajectory (TIK+TV)",'Interpreter','Latex')

title("Gaussianity 3d plot")


figure

[mu2_ca,mu1_ca]=find(W==min(W(:)));

contourf(log(MU1),log(MU2),log(W),'LineWidth',1,'EdgeColor', 'none')
hold on
scatter(log(mu1(mu1_ca)),log(mu2(mu2_ca)),'r o','LineWidth',1)
hold on
scatter(betaa(k(1)),betaa(l(1)),'+','LineWidth',1)
hold on
scatter(betaa2(k(1)),betaa2(l(1)),'+','LineWidth',1)
hold on
plot(betaavec(1,:),betaavec(2,:),':','LineWidth',1)
hold on
plot(betaavec2(1,:),betaavec2(2,:),':','LineWidth',1)


legend("$W=\widetilde{W}(\mu_1$,$\mu_2)$","$W=\widetilde{W}(\mu_1^*,\mu_2^*)$ (grid)","$W=\widetilde{W}(\mu_1^*,\mu_2^*)$ (GN)","$W=\widetilde{W}(\mu_1^*,\mu_2^*)$ (GN TIK+TV)","GN trajectory","GN trajectory (GN TIK+TV)",'Interpreter','Latex')

xlabel("log($\mu_1$)",'Interpreter','Latex')
ylabel("log($\mu_2$)",'Interpreter','Latex')

title("Whiteness level curve")

figure

surfc(log(MU1),log(MU2),log(W),'LineWidth',1,'EdgeColor', 'none')
%hold on
%contour3(log(MU1),log(MU2),log(g),'LineWidth',1)
hold on
plot3(log(mu1(mu1_ca)),log(mu2(mu2_ca)),log(W(mu2_ca,mu1_ca)),'r o','LineWidth',1)
hold on
plot3((betaa(k(1))),(betaa(l(1))),log(fk3(end)),'*','LineWidth',1)
hold on
plot3((betaa2(k(1))),(betaa2(l(1))),log(fk32(end)),'g o','LineWidth',1)
hold on
plot3((betaavec(1,:)),(betaavec(2,:)),log(fk3),'k .- ','LineWidth',1)
hold on
plot3((betaavec2(1,:)),(betaavec2(2,:)),log(fk32),'b .- ','LineWidth',1)

legend("$W=\widetilde{W}(\mu_1,\mu_2)$","$W=\widetilde{W}(\mu_1,\mu_2)$ level curve","$W=\widetilde{W}(\mu_1^*,\mu_2^*)$ (grid)","$W=\widetilde{W}(\mu_1^*,\mu_2^*)$ (GN)","$W=\widetilde{W}(\mu_1^*,\mu_2^*)$ (GN TIK + TV)","GN trajectory","GN trajectory (TIK+TV)",'Interpreter','Latex')

[q2,q1]=find(PSNR==max(PSNR(:)));

figure
contourf(log(MU1),log(MU2),PSNR,'Edgecolor','none')
hold on
scatter(log(mu1(q1)),log(mu2(q2)),'+','LineWidth',1)
hold on
scatter(betaa_err(k(1)),betaa_err(l(1)),'Linewidth',1)
hold on
scatter(betaa_(k(1)),betaa_(l(1)),'Linewidth',1)
hold on
scatter(betaa(k(1)),betaa(l(1)),'Linewidth',1)
legend("PSNR level curve",'max PSNR','MSE','gaussianity','Whiteness')

[q2,q1]=find(SSIM==max(SSIM(:)));


figure
contourf(log(MU1),log(MU2),SSIM,'Edgecolor','none')
hold on
scatter(log(mu1(q1)),log(mu2(q2)),'+','LineWidth',1)
hold on
scatter(betaa_err(k(1)),betaa_err(l(1)),'Linewidth',1)
hold on
scatter(betaa_(k(1)),betaa_(l(1)),'Linewidth',1)
hold on
scatter(betaa(k(1)),betaa(l(1)),'Linewidth',1)
legend("PSNR level curve",'max SSIM','MSE','gaussianity','Whiteness')

figure
subplot(1,3,1)
imshow2(xF)
subplot(1,3,2)
imshow2(bb)
subplot(1,3,3)
imshow2(M)

figure
subplot(1,3,1)
imshow2(x1)
title("Deblurred image with MSE")
subplot(1,3,2)
imshow2(x2)
title("Deblurred image with Gaussianity")
subplot(1,3,3)
imshow2(x3)
title("Deblurred image with Whiteness")






