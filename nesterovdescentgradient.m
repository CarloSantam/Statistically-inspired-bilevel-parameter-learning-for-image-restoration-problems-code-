function [x,ff]=nesterovdescentgradient(maxit,x0,mu,bb,epsi,H_FT,tol)
x=x0;
x1=x;
%y=y0;
t=1;

alpha=1/((mu)+12/epsi);
[f2]=fun(mu,bb,epsi,x,H_FT);
ff(1)=f2;

for i=1:maxit
    
    tnew=(1+sqrt(1+4*t^2))/2;
    y=x+((t-1)/tnew)*(x-x1);

    [g]=gradient(H_FT,y,bb,epsi,mu);
    z=y-alpha*g;

    x1=x;

    x=z;

    

    t=tnew;

    [f2]=fun(mu,bb,epsi,x,H_FT);

    ff(i+1)=f2;

    if norm(ff(i+1)-ff(i))/norm(ff(i))<tol
        break
    end

end
end

function [f]=fun(mu,bb,epsi,x,H_FT)

z=real(ifft2(H_FT.*fft2(x)))-bb;
Dhx=[ x(:,2:end) - x(:,1:(end-1)) , x(:,1) - x(:,end) ];
Dvx=[ x(2:end,:) - x(1:(end-1),:) ; x(1,:) - x(end,:) ];

if 0
[L1,dL1]=huber(Dhx(:),epsi);


[L2,dL2]=huber(Dvx(:),epsi);

f=mu/2*norm(z,'fro')^2+sum(L1)+sum(L2);

dL1=reshape(dL1,[length(Dhx),length(Dhx)]);
dL2=reshape(dL2,[length(Dvx),length(Dvx)]);

else

    [L1]=huber_0(Dhx,epsi);

    [L2]=huber_0(Dvx,epsi);

    

    f=(mu/2)*norm(z,'fro')^2+sum(L1(:))+sum(L2(:));

end

end

function g=gradient(H_FT,x,bb,epsi,mu)

z=real(ifft2(H_FT.*fft2(x)))-bb;

Dhx=[ x(:,2:end) - x(:,1:(end-1)) , x(:,1) - x(:,end) ];
Dvx=[ x(2:end,:) - x(1:(end-1),:) ; x(1,:) - x(end,:) ];
%DhTx = [ x(:,end) - x(:,1) , -diff(x,1,2) ];
%DvTx = [ x(end,:) - x(1,:) ; -diff(x,1,1) ];

[dL1]=huber_1(Dhx,epsi);

[dL2]=huber_1(Dvx,epsi);

g1=mu*real(ifft2(conj(H_FT).*fft2(z))); %mu*H^T(Hx-b)

g2=[ dL1(:,end) - dL1(:,1) , -diff(dL1,1,2) ];

g3=[ dL2(end,:) - dL2(1,:) ; -diff(dL2,1,1) ];


g=g1+g2+g3;

end

