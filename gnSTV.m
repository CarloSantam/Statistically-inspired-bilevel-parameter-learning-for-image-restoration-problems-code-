function [betaa,fk]=gnSTV(type,Dh_FT,Dv_FT,x0,z0,H_FT,bb,xF,...
    beta_0,maxit1,maxit2,tol,alpha,epsi,alphaa,beta,sigma,n)

betaa=beta_0;

for i=1:maxit1
    [Jac,f]=nonsmoothparametergradfun(type,maxit2,x0,z0,alphaa,beta,betaa,bb,Dh_FT,Dv_FT,H_FT,tol,xF,epsi,sigma,n);
    d=-(Jac'*f)/(Jac'*Jac)
    if norm(d)<tol
        break
    end
    betaa=betaa+alpha*d
    [Jac,f]=nonsmoothparametergradfun(type,maxit2,x0,z0,alphaa,beta,betaa,bb,Dh_FT,Dv_FT,H_FT,tol,xF,epsi,sigma,n);
    fk(i)=norm(f)^2;
    
end
end