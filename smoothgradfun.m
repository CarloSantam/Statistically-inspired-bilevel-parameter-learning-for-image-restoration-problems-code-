function [f,g]=smoothgradfun(mu,bb,x,Dh_FT,Dv_FT,H_FT)

z=real(ifft2(H_FT.*fft2(x)))-bb;

Dhx=real(ifft2(Dh_FT.*fft2(x)));
Dvx=real(ifft2(Dv_FT.*fft2(x)));

f=mu*norm(z,'fro')^2+norm(Dhx,'fro')^2+norm(Dvx,'fro')^2;
g1=2*mu*real(ifft2(conj(H_FT).*fft2(z)));
g2=2*real(ifft2(conj(Dh_FT).*fft2(Dhx)));
g3=2*real(ifft2(conj(Dv_FT).*fft2(Dvx)));
g=g1+g2+g3;
end