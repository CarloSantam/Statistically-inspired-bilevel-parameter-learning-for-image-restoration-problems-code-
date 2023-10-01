function [beta,fk,time]=gn1D(DTD_FT,H_FT,HTH_FT,bbhat,xF,beta_0,maxit,tol,type,sigma,n,bb)

beta=beta_0;

time(1)=0;

if type==1

for i=1:maxit
    tstart=tic;
    [grad,f]=gradfun(DTD_FT,H_FT,HTH_FT,beta,bbhat,xF);
    d=-(grad'*f)/(grad'*grad);
    fk(i)=norm(f)^2;
    if norm(d)<tol
        break
    end
    beta=beta+d;
    a=toc(tstart);
    time(i+1)=time(i)+a;
end

elseif type==2

for i=1:maxit
    tstart=tic;
    [grad,f]=gradfungauss(DTD_FT,H_FT,HTH_FT,beta,bbhat,n,sigma);
    d=-(grad'*f)/(grad'*grad);
    fk(i)=norm(f)^2;
    if norm(d)<tol
        break
    end
    beta=beta+d;
    a=toc(tstart);
    time(i+1)=a+time(i);
end

elseif type==3

for i=1:maxit
    tstart=tic;
    [grad,f]=gradfunwhiteness(DTD_FT,H_FT,HTH_FT,beta,bbhat,bb);
    d=-(grad'*f)/(grad'*grad);
    fk(i)=norm(f)^2;
    if norm(d)<tol
        break
    end
    beta=beta+d;
    a=toc(tstart);
    time(i+1)=a+time(i);
end

else
    disp("This type is not avaiable")
end
end
