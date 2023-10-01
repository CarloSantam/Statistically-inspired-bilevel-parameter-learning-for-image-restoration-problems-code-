function [x,ff,kk,alphavec,ymatrix]=smoothnesterovdescentgradient(maxit,x0,y0,alpha_0,mu,bb,Dh_FT,Dv_FT,H_FT,c,tau,tol,xFF)
x=x0;
y=y0;
t=1;
ymatrix=[];

for i=1:maxit
    x1=x;
    alpha=alpha_0;
    [f1,g]=smoothgradfun(mu,bb,y,Dh_FT,Dv_FT,H_FT);
    m=-g(:)'*g(:);
    tt=-c*m;
    x_new=y-alpha*g;
    [f2,~]=smoothgradfun(mu,bb,x_new,Dh_FT,Dv_FT,H_FT);
    while f1-f2<alpha*tt
        alpha=alpha*tau;
        x_new=y-alpha*g;
        [f2,~]=smoothgradfun(mu,bb,x_new,Dh_FT,Dv_FT,H_FT);    
    end
    ff(i)=f2;
    x=x_new;
    if norm(x1(:)-x(:))<tol
        break
    end
    tnew=(1+sqrt(1+4*t^2))/2;
    y=x+((t-1)/tnew)*(x-x1);
    t=tnew;
    kk(i)=norm(x-xFF,'fro');
    ymatrix=[ymatrix,y(:)];
    alphavec(i)=alpha;
end