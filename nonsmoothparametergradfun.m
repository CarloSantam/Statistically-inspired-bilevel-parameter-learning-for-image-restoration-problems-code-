function [grad,fun]=nonsmoothparametergradfun(type,maxit,x0,z_0,alpha,beta,betaa,...
    bb,Dh_FT,Dv_FT,H_FT,tol,xF,epsi,sigma,n)
mu=exp(betaa);
[x,~,ymatrix,~,~ ,iter]=momentumdescentgradient(maxit,x0,z_0,alpha,beta,mu,bb,epsi,Dh_FT,Dv_FT,H_FT,tol);

ck=length(H_FT);
s=zeros(ck,ck);
for i=0:(iter-1)
    z=real(ifft2(H_FT.*fft2(reshape(ymatrix(:,i+1),[ck,ck]))))-bb;
    l=real(ifft2(conj(H_FT).*fft2(z)));
    s=s-2*mu*alpha*((1-beta^(iter-i))/(1-beta))*l;
end
if type==1
    grad=s(:);
    fun=x(:)-xF(:);
elseif type==2
    grad=real(ifft2(H_FT.*fft2(s)));
    Res=real(ifft2(H_FT.*fft2(x)))-bb;
    grad=grad(:);
    Res=Res(:);
    grad=2*grad'*Res;
    fun=norm(Res)^2-(sigma^2)*n^2;
elseif type==3
    gder=real(ifft2(H_FT.*fft2(s)));
    Res=real(ifft2(H_FT.*fft2(x)))-bb;
    ccorrelation1=real(ifft2(fft2(gder).*conj(fft2(Res))))+real(ifft2(fft2(Res).*conj(fft2(gder))));
    ccorrelation2=real(ifft2(fft2(Res).*conj(fft2(Res))));
    gder=gder(:);
    Res=Res(:);
    Resnorm=norm(Res,2);
    resgrad=2*gder'*Res;
    fun=ccorrelation2(:)/Resnorm^2;
    grad=(ccorrelation1(:)*Resnorm^2-ccorrelation2(:)*resgrad)/Resnorm^4;
else
    disp("This thype is not avaiable")
end

end