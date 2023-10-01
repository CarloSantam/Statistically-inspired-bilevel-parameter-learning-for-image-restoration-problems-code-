function [grad,fun,ss,relres]=nonsmoothgradientimplicit(type,epsi,H_FT,HTH_FT,beta,...
    xstar,xF,bb,maxit,tol,sigma,ss0)

mu=exp(beta);

muHTH_FT=mu*HTH_FT;

%minustwomuHTb=mu*minustwoHTb;

RHS=(-mu)*real(ifft2 ( conj(H_FT).*(H_FT.*fft2(xstar) - fft2(bb) ) ) );

RHS=RHS(:);

Dhu=Dh(xstar);
Dvu=Dv(xstar);

if 0
%[~,~,dd1]=huber(Dhu(:),epsi);
%[~,~,dd2]=huber(Dvu(:),epsi);
%dd1=reshape(dd1,[length(Dhu),length(Dhu)]);
%dd2=reshape(dd2,[length(Dvu),length(Dvu)]);
else

    [dd1]=huber_2(Dhu,epsi);
    [dd2]=huber_2(Dvu,epsi);
end

n=length(H_FT);


[ss,~,relres]=minres(@(x)hessian(x,muHTH_FT,dd1,dd2),RHS,tol,maxit,[],[],ss0);

s=reshape(ss,[n,n]);

funn=xstar(:)-xF(:);
if type==1
    grad=s(:);
    fun=funn;
elseif type==2
    grad=real(ifft2(H_FT.*fft2(s)));
    Res=real(ifft2(H_FT.*fft2(xstar)))-bb;
    grad=grad(:);
    Res=Res(:);
    grad=2*grad'*Res;
    fun=norm(Res)^2-(sigma^2)*n^2;
elseif type==3
    gder=real(ifft2(H_FT.*fft2(s)));
    Res=real(ifft2(H_FT.*fft2(xstar)))-bb;
    ccorrelation1=real(ifft2(fft2(gder).*conj(fft2(Res))))+real(ifft2(fft2(Res).*conj(fft2(gder))));
    ccorrelation2=real(ifft2(fft2(Res).*conj(fft2(Res))));
    gder=gder(:);
    Res=Res(:);
    Resnorm=norm(Res,2);
    resgrad=2*gder'*Res;
    fun=ccorrelation2(:)/Resnorm^2;
    grad=(ccorrelation1(:)*Resnorm^2-ccorrelation2(:)*resgrad)/Resnorm^4;
end
end

function HE = hessian(x,muHTH_FT,dd1,dd2)
x=reshape(x,[size(dd1)]);
HE=real(ifft2(muHTH_FT.*fft2(x)))+DhT(dd1.*Dh(x))+DvT(dd2.*Dv(x));
HE=HE(:);
end

function v = vec(u)
  v = u(:);
end
function Dhu = Dh(u)
    Dhu  = [ u(:,2:end) - u(:,1:(end-1)) , u(:,1) - u(:,end) ];
end
function Dvu = Dv(u)
    Dvu  = [ u(2:end,:) - u(1:(end-1),:) ; u(1,:) - u(end,:) ]; 
end
function DhTu = DhT(u)
    DhTu = [ u(:,end) - u(:,1) , -diff(u,1,2) ];
end
function DvTu = DvT(u)
    DvTu = [ u(end,:) - u(1,:) ; -diff(u,1,1) ];
end
