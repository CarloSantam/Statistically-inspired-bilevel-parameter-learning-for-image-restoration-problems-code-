function [x,ff,xmatrix,alpha,beta,iter,xx]=momentumdescentgradient(maxit,x_0,z_0,alpha,beta,mu,bb,epsi,Dh_FT,Dv_FT,H_FT,tol)
x=x_0;
xmatrix=[x_0(:)];
z=z_0;
k=0;
for i=1:maxit
    [f1,g]=nonsmoothgradfun(mu,bb,epsi,x,Dh_FT,Dv_FT,H_FT);
    
    if norm(g(:))<tol
        break
    end
    
    z=beta*z+g;
    x=x-alpha*z;
    
    k=k+1;
    xmatrix=[xmatrix,x(:)];
    [f1,g]=nonsmoothgradfun(mu,bb,epsi,x,Dh_FT,Dv_FT,H_FT);
    ff(i)=f1;
end
iter=k;

c=length(H_FT);
xx=x_0;

for i=0:(iter-1)
    [~,g]=nonsmoothgradfun(mu,bb,epsi,reshape(xmatrix(:,i+1),[c,c]),Dh_FT,Dv_FT,H_FT);
    kk=alpha*((1-beta^(iter-i))/(1-beta))*g;
    xx=xx-kk;
end

end
