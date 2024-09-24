function [x,tau] = VAMP_MMV_2d_SVD_c(y_in,H,N0,prow,pcol,vpath)%-----------
    [Nr, Nt] = size(H);
    [~,Nall] = size(y_in);
    x_hat = zeros(Nt,Nall); %(x(-1))
    r_t = zeros(Nt,Nall);
    gamma = prow*ones(1,Nall);
    gamma_t = Nr/Nt*ones(1,Nall);
    % gamma_t = 1*ones(1,Nall);
    v_hat = ones(1,Nall);

    [U,S,V]=svd(H);
    s = diag(S);
    R = min(size(H));
    V = V(:,1:R);
    U = U(:,1:R);
    y_t = (diag(s)\(U'))*y_in;

    damping = 0.8;
    tol = 2e-6;

    % Eg2 = @(A,y,r2,N,ga2)double((A'*A/N+ga2*eye(size(A,2)))\(A'*y/N+ga2*r2));
    % Vg2 = @(A,y,r2,N,ga2)double(trace((A'*A/N+ga2*eye(size(A,2)))^(-1)))/(Nt);
    dfunc = @(s,N0,ga) ((1/N0)./(1./N0*abs(s).^2+ga)).*abs(s).^2;
    Iternum = 18;
    for k = 1:Iternum
        for kant = 1:Nall  
            d(:,kant) = dfunc(s,N0,gamma_t(kant));
            gamma(kant) = gamma_t(kant)*mean(d(:,kant))/(Nt/R-mean(d(:,kant)));
            gamma(kant) = min(gamma(kant),1e6);
            r(:,kant) = r_t(:,kant) + Nt/R*V*diag(d(:,kant)/mean(d(:,kant)))*(y_t(:,kant)-V'*r_t(:,kant));
        end

        v_hat_pre = v_hat;
        [x_hat_temp,v_hat] = g1_MMV_2d(r,gamma,prow,pcol,vpath);
        for kant = 1:Nall
            alpha(kant) = v_hat(kant)*gamma(kant);
            gamma_t_pre(kant) = gamma_t(kant);
            gamma_t(kant) = gamma(kant)*(1-alpha(kant))/alpha(kant);
            gamma_t(kant) = min(gamma_t(kant),1e5);
            if gamma_t(kant)<0
                gamma_t(kant)=gamma_t_pre(kant);
                x_hat(:,kant) = x_hat(:,kant);
            else 
                x_hat(:,kant) = x_hat_temp(:,kant)*damping+(1-damping)*x_hat(:,kant);
            end
            r_t(:,kant) = (x_hat(:,kant)-alpha(kant)*r(:,kant))/(1-alpha(:,kant));
        end
        if max(abs(v_hat_pre-v_hat))<1e-3
            break;
        end
    end
    x = x_hat;
    tau = v_hat;
end


function [miu,sigma] = g1_MMV_2d(r,gamma,prow,pcol,vpath)
    [Nt,ndim] = size(r);
    CN = @(x,r,sigma)(1./(pi^ndim*prod(sigma)).*exp(-(x-r)'*diag(1./sigma)*(x-r)));
    logCN = @(x,r,sigma)(-(x-r)'*diag(1./sigma)*(x-r)-log((pi^ndim*prod(sigma))));
    CN1 = @(x,r,sigma)(1./(pi*sigma).*exp(-conj(x-r)*diag(1/sigma)*(x-r)));
    sigma_priori = vpath;
    for k1 = 1:Nt
        for k2 = 1:ndim
            pnum(k1,k2) =(1-pcol)*real(CN1(0,r(k1,k2),1/gamma(k2)));
            pden(k1,k2) = pcol*real(CN1(0,r(k1,k2)-0,sigma_priori+1/gamma(k2)));
            % pnum(k1,k2) = max(pnum(k1,k2),1e-4);
            % pden(k1,k2) = max(pden(k1,k2),1e-4);
            pval(k1,k2) = (1-pcol)/pcol*(sigma_priori+1/gamma(k2))*gamma(k2)*exp(-abs(r(k1,k2)).^2*(gamma(k2)-1/(sigma_priori+1/gamma(k2))));
            sigma_out(k1,k2) = 1/(1/sigma_priori+gamma(k2));
            rout(k1,k2) = sigma_out(k1,k2)*(r(k1,k2)*gamma(k2)); 
        end
        normdata = min(CN(0,r(k1,:).',1./gamma)/prod((pnum(k1,:)+pden(k1,:)))*(1-prow),1e300)+prow;
        for k2 = 1:ndim
            miu(k1,k2) = prow/(1+pval(k1,k2))*rout(k1,k2)/normdata;
            sigma_all(k1,k2) = prow/(1+pval(k1,k2))*(sigma_out(k1,k2)+abs(rout(k1,k2)).^2)/normdata-abs(miu(k1,k2)).^2;
        end
    end
    sigma = mean(abs(sigma_all));
end