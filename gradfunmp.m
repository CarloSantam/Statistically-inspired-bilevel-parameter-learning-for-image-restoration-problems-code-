function [f,g]=gradfunmp(H_FT,x,bb,mu,epsi)

z=real(ifft2(H_FT.*fft2(x)))-bb;
Dhx=Dh(x);
Dvx=Dv(x);

[n,~]=size(bb);
mu_mat = ones(size(bb));
mu_mat(1:n,(n/2+1):n)=mu(2);
mu_mat(1:n,1:n/2)=mu(1);

[L1,dL1]=huber(Dhx(:),epsi);

[L2,dL2]=huber(Dvx(:),epsi);

A=abs(z).^2;

f=sum((mu_mat(:)).*A(:))+sum(L1+L2);

dL1=reshape(dL1',[size(Dhx)]);
dL2=reshape(dL2',[size(Dhx)]);

%s=real(ifft2(conj(H_FT).*fft2(z)));
s=(2*mu_mat).*z;
g1=real(ifft2(conj(H_FT).*fft2(s)));
 %2*mu.*H^T(Hx-b)
g2=DhT(dL1);
g3=DvT(dL2);

g=g1+g2+g3;
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
