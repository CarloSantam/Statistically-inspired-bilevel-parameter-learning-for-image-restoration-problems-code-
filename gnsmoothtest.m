function [beta,fk]=gnsmoothtest(Dh_FT,Dv_FT,x0,y0,H_FT,bb,xF,xFF,...
    beta_0,maxit1,maxit2,tol,tol1,c,tau,c1,tau1,alpha_0)

beta=beta_0;

for i=1:maxit1
    alpha=alpha_0;
    [Jac,f,Sgrad]=smoothparametergradfun(maxit2,x0,y0,alpha_0,beta,bb,Dh_FT,Dv_FT,H_FT,c,tau,tol,xFF,xF);
    d=-(Jac'*f)/(Jac'*Jac);
    m=Sgrad'*d;
    if norm(d)<tol
        break
    end
    t=-c1*m;
    beta_new=beta+alpha*d;
    f1=norm(f)^2;
    [~,f,~]=smoothparametergradfun(maxit2,x0,y0,alpha_0,beta_new,bb,Dh_FT,Dv_FT,H_FT,c,tau,tol,xFF,xF);
    f2=norm(f)^2;
    while f1-f2<alpha*t
        alpha=tau1*alpha
        if alpha<tol1
            break
        end
        beta_new=beta+alpha*d;
        [~,f,~]=smoothparametergradfun(maxit2,x0,y0,alpha_0,beta_new,bb,Dh_FT,Dv_FT,H_FT,c,tau,tol,xFF,xF);
        f2=norm(f)^2;
    end
    
    
    beta=beta_new
    fk(i)=f2;
end


end