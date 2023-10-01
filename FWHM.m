clear all
close all
clc

sigma=3;
support=41;
[gaussian,center]=psfGauss(support,sigma);

% Calcola la FWHM (Full Width at Half Maximum) basata sulle deviazioni standard
fwhm_x = 2.355 * sigma;
fwhm_y = 2.355 * sigma;

% Calcola il valore della gaussiana 2D
%x = linspace(x_mean - 2*fwhm_x, x_mean + 2*fwhm_x, 100);
%y = linspace(y_mean - 2*fwhm_y, y_mean + 2*fwhm_y, 100);
%[X, Y] = meshgrid(x, y);
%gaussian = exp(-((X - x_mean).^2 / (2 * sigma^2) + (Y - y_mean).^2 / (2 * sigma^2)));

% Calcola l'altezza a cui vuoi posizionare il quadrato (metà dell'altezza massima)
half_height = max(gaussian(:)) / 2;

% Calcola la larghezza del quadrato
width = fwhm_x;
height = fwhm_y;

% Calcola le coordinate del rettangolo
%x_rect = center(1) - width / 2;
%y_rect = center(2) - height / 2;

% Crea una figura vuota
figure;

% Visualizza la gaussiana 2D
surf(gaussian);
xlabel('x');
ylabel('y');
title('Gaussian PSF with $FWHM_x$ and $FWHM_y$ in evidence','Interpreter','latex');

[K,J]=meshgrid(center(1) - (width / 2):0.01:center(1) + (width / 2),...
    center(2) - (height / 2):0.01:center(2) + (height / 2));
hold on
% Disegna il rettangolo
plot3(K(1,:),J(1,:),half_height*ones(length(K),1),'b','Linewidth',2)
hold on
plot3(K(:,1),J(:,1),half_height*ones(length(K),1),'k','Linewidth',2)
hold on
plot3(K(end,:),J(end,:),half_height*ones(length(K),1),'b','Linewidth',2)
hold on
plot3(K(:,end),J(:,end),half_height*ones(length(K),1),'k','Linewidth',2)
legend("$PSF$","$FWHM_x$","$FWHM_y$","","","interpreter","latex")

figure

n=41;
y = -fix(n/2):0.00001:ceil(n/2)-1;
for i=1:length(y)
a(i)=exp(-(y(i).^2)./(2*sigma^2));
end

aa=a/sum(a(:));

plot(y,aa,'LineWidth',1)
axis tight;
hold on
max_h=max(aa(:))/2
x=max_h;

a=-fwhm_x/2:0.0000001:fwhm_x/2;

plot(a,max_h*ones(size(a)),'r','LineWidth',1)

ylabel("$y$",'interpreter','latex')
xlabel("$x$",'interpreter','latex')

legend("PSF 1D", "FWHM",'interpreter','latex')
