
function [betaa,fk,xstar,time,PNSR,SSIM]=gnSTVimplicit(type,maxit,maxit1,x_0,alpha,beta_0,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,sigma,tol_r,toln,rel_res_sf)
    
betaa=beta_0;

HTH_FT=(conj(H_FT).*H_FT);

%minustwoHTb=-2*real(ifft2(conj(H_FT).*fft2(bb)));

%[Jac,f,ss]=nonsmoothgradientimplicit(type,epsi,H_FT,HTH_FT,betaa,...
    %x_0,xF,bb,maxitpcg,tolpcg,sigma,[]);

fk=nan(maxit1,1);

time=nan(maxit1,1);

time(1)=0;

PNSR=nan(maxit1,1);

SSIM=nan(maxit1,1);

xstar=x_0;

ss=[];

for i=1:maxit1
 
    tstart=tic;
    xstar=nesterovdescentgradient(maxit,xstar,exp(betaa),bb,epsi,H_FT,toln);
    [Jac,f,ss,relres]=nonsmoothgradientimplicit(type,epsi,H_FT,HTH_FT,betaa,xstar,xF,bb,maxitpcg,tolpcg,sigma,ss);
    
    normf2=norm(f)^2;

    tolpcg=relres/rel_res_sf;
    
    d=(Jac'*f)/(-Jac'*Jac);
    
    betaa_new=betaa+alpha*d

    a=toc(tstart);
    
    time(i+1)=time(i)+a;
    fk(i+1)=norm(f)^2;
    PNSR(i+1)=psnr(xF,xstar);
    SSIM(i+1)=ssim(xF,xstar);


    if norm(normf2-fk(i))/norm(normf2)<tol_r || norm(d)<tol
        break
    end

    betaa=betaa_new;
    
end

time=rmmissing(time(2:end));
fk=rmmissing(fk);
PNSR=rmmissing(PNSR);
SSIM=rmmissing(SSIM);

end
