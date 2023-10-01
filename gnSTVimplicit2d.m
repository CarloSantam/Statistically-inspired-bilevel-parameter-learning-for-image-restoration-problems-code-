function [betaa,fk,xstar,time,PNSR,SSIM,betaavec,ff]=gnSTVimplicit2d(type,maxit,maxit1,x0,alpha,beta_0,bb,...
    epsi,H_FT,tol,xF,maxitpcg,tolpcg,tol_r,toln,sigma1,sigma2,sigma_mat,rel_res_sf,M,k,l)

    time=[];

    %twoHTH_FT=2*(conj(H_FT).*H_FT);

    %betaa=ones(n);

    %betaa(1:n,1:n/2)=beta_0(1);

    %betaa(1:n,(n/2+1):end)=beta_0(2);

    betaa=beta_0(1)*M+beta_0(2)*(ones(size(M))-M);    
    
    %[Jac,f,ss1,ss2]=nonsmoothgradientimplicit2p(type,epsi,H_FT,betaa,x0,xF,bb,maxitpcg,tolpcg,[],[],sigma1,sigma2,M,sigma_mat,k,l);
    
    fk=nan(maxit1,1);
    time=nan(maxit1,1);
    time(1)=0;
    PNSR=nan(maxit1,1);
    SSIM=nan(maxit1,1);

    xstar=x0;

    ss1=[];
    ss2=[];

    betaavec=nan(2,maxit1);
    

    for i=1:maxit1


        tstart=tic;

        [xstar,ff]=nesterovdescentgradientmp(maxit,xstar,exp(betaa),bb,epsi,H_FT,toln);
       
        [Jac,f,ss1,ss2,resf,resf2]=nonsmoothgradientimplicit2p(type,epsi,H_FT,betaa,xstar,xF,bb,maxitpcg,tolpcg,ss1,ss2,sigma1,sigma2,M,sigma_mat,k,l);
        betaavec(:,i+1)=[betaa(k(1));betaa(l(1))];

        normf2=norm(f)^2
        
        d=(Jac'*Jac)\(-Jac'*f)

        dd=d(1)*M+d(2)*(ones(size(M))-M);

        betaa_new=betaa+alpha*dd;

         %if norm(normf2-fk(i))/norm(normf2)<tol_r || norm(d)<tol
         %   break
         %end
         
        %betaavec(:,i+1)=[betaa(k(1));betaa(l(1))];

        tolpcg=min(resf/rel_res_sf,resf2/rel_res_sf);
                
        %dd=d(1)*M+d(2)*(ones(size(M))-M);

        %betaa=betaa+alpha*dd;
        
        a=toc(tstart);
        
        time(i+1)=time(i)+a;
        PNSR(i+1)=psnr(xF,xstar);
        SSIM(i+1)=ssim(xF,xstar);
        fk(i+1)=normf2;

        if norm(normf2-fk(i))/norm(normf2)<tol_r || norm(d)<tol
            break
        end

        betaa=betaa_new;
        
    end

    time=rmmissing(time(2:end));
    PNSR=rmmissing(PNSR);
    SSIM=rmmissing(SSIM);
    fk=rmmissing(fk);
    betaavec=rmmissing(betaavec,2);
    
end
