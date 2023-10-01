function [grad,fun,Sgrad]=smoothparametergradfun(maxit,x0,y0,alpha_0,beta,bb,Dh_FT,Dv_FT,H_FT,c,tau,tol,xFF,xF)
mu=exp(beta);
[x,~,~,alphavec,ymatrix]=smoothnesterovdescentgradient(maxit,x0,y0,alpha_0,mu,bb,Dh_FT,Dv_FT,H_FT,c,tau,tol,xFF);

ck=length(H_FT);
s=zeros(length(H_FT)^2,1);
for i=1:length(alphavec)
    z=real(ifft2(H_FT.*fft2(reshape(ymatrix(:,i),[ck,ck]))))-bb;
    l=real(ifft2(conj(H_FT).*fft2(z)));
    k=alphavec(i)*l(:);
    s=s+k;

end
grad=-mu*s;
fun=x(:)-xF(:);
Sgrad=2*grad'*fun;
end