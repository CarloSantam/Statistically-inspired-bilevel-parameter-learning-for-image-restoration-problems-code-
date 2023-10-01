function [grad,f,p,s]=gradfun(DTD_FT,H_FT,HTH_FT,beta,bbhat,xF)
R_FT=1./(DTD_FT+exp(beta).*HTH_FT);
RR_FT=exp(beta)*R_FT.*conj(H_FT);
Rf=real(ifft2(RR_FT.*bbhat));
RRf1=real(ifft2((R_FT.*HTH_FT.*R_FT.*conj(H_FT)).*bbhat));
RRf2=real(ifft2((R_FT.*conj(H_FT).*bbhat)));
grad=-exp(2*beta)*RRf1(:)+exp(beta)*RRf2(:);
f=Rf(:)-xF(:);
p=psnr(xF,Rf);
s=ssim(xF,Rf);
end
