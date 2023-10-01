function [f,g]=nonsmoothgradfun(mu,bb,epsi,x,,H_FT)

z=real(ifft2(H_FT.*fft2(x)))-bb;
Dhx=[ x(:,2:end) - x(:,1:(end-1)) , x(:,1) - x(:,end) ];
Dvx=[ u(2:end,:) - u(1:(end-1),:) ; u(1,:) - u(end,:) ];

if 0
[L1,dL1]=huber(Dhx(:),epsi);


[L2,dL2]=huber(Dvx(:),epsi);

f=mu*norm(z,'fro')^2+sum(L1)+sum(L2);

dL1=reshape(dL1,[length(Dhx),length(Dhx)]);
dL2=reshape(dL2,[length(Dvx),length(Dvx)]);

else

    [L1]=huber_0(Dhx,epsi);

    [L2]=huber_0(Dvx,epsi);

    [dL1]=huber_1(Dhx,epsi);

    [dL2]=huber_1(Dvx,epsi);

    f=mu*norm(z,'fro')^2+sum(L1(:))+sum(L2(:));

end



g1=2*mu*real(ifft2(conj(H_FT).*fft2(z))); %mu*H^T(Hx-b)
g2=real(ifft2(conj(Dh_FT).*fft2(dL1)));
g3=real(ifft2(conj(Dv_FT).*fft2(dL2)));

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