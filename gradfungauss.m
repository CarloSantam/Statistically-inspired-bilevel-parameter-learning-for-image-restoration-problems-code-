function [grad,f]=gradfungauss(DTD_FT,H_FT,HTH_FT,beta,bbhat,n,sigma)
R_FT=1./(DTD_FT+exp(beta).*HTH_FT);
RR_FT=exp(beta)*R_FT.*conj(H_FT);
Rf=real(ifft2(RR_FT.*bbhat));
RRf1=real(ifft2((R_FT.*HTH_FT.*R_FT.*conj(H_FT)).*bbhat));
RRf2=real(ifft2((R_FT.*conj(H_FT).*bbhat)));
grad=-exp(2*beta)*RRf1+exp(beta)*RRf2;

gder=real(ifft2(H_FT.*fft2(grad)));
Res=real(ifft2(H_FT.*fft2(Rf)))-ifft2(bbhat);
gder=gder(:);
Res=Res(:);

grad=2*gder'*Res;
f=norm(Res)^2-(sigma^2)*n^2;

end