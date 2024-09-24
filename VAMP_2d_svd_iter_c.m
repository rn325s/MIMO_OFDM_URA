function [x,tau] = VAMP_2d_svd_iter_c(y_in,H_all,N0,ModOrder,H_LDPC,rndseed)%----------------
    [Nr, Nt,Nslot] = size(H_all);
    for k = 1:Nslot
        H{k} = reshape(H_all(:,:,k),Nr,Nt);
        [U1,S1,V1]=svd(H{k});
        if Nt == 1
            S1 = S1(1);
        end
        s{k} = diag(S1);
        R = min(size(H{k}));
        V{k} = V1(:,1:R);
        U{k} = U1(:,1:R);
        y_t(:,k) = (diag(s{k})\(U{k}'))*y_in(:,k);
    end
    x_hat = zeros(Nt,Nslot); %(x(-1))
    r_t = zeros(Nt,Nslot);
    gamma = ones(1,Nslot);
    gamma_t = 1*ones(1,Nslot);



    damping = 0.7;
    tol = 2e-6;

    % Eg2 = @(A,y,r2,N,ga2)double((A'*A/N+ga2*eye(size(A,2)))\(A'*y/N+ga2*r2));
    % Vg2 = @(A,y,r2,N,ga2)double(trace((A'*A/N+ga2*eye(size(A,2)))^(-1)))/(Nt);
    dfunc = @(s,N0,ga) ((1/N0)./(1./N0*abs(s).^2+ga)).*abs(s).^2;
    Iternum = 18;
    for k = 1:Iternum
        for kslot = 1:Nslot  
            d(:,kslot) = dfunc(s{kslot},N0,gamma_t(kslot));
            gamma(kslot) = gamma_t(kslot)*mean(d(:,kslot))/(Nt/R-mean(d(:,kslot)));
            gamma(kslot) = min(gamma(kslot),1e5);
            r(:,kslot) = r_t(:,kslot) + Nt/R*V{kslot}*diag(d(:,kslot)/mean(d(:,kslot)))*(y_t(:,kslot)-V{kslot}'*r_t(:,kslot));
        end

        [x_hat_temp,v_hat] = g1_iter_2d_mu(r,gamma,ModOrder,H_LDPC,rndseed,Nt);
        for kslot = 1:Nslot
            alpha(kslot) = v_hat(kslot)*gamma(kslot);
            gamma_t_pre(kslot) = gamma_t(kslot);
            gamma_t(kslot) = gamma(kslot)*(1-alpha(kslot))/alpha(kslot);
            gamma_t(kslot) = min(gamma_t(kslot),1e5);
            if gamma_t(kslot)<0
                gamma_t(kslot)=gamma_t_pre(kslot);
                x_hat(:,kslot) = x_hat(:,kslot);
            else 
                x_hat(:,kslot) = x_hat_temp(:,kslot)*damping+(1-damping)*x_hat(:,kslot);
            end
            r_t(:,kslot) = (x_hat(:,kslot)-alpha(kslot)*r(:,kslot))/(1-alpha(:,kslot));
        end
        
    end
    x = x_hat;
    tau = v_hat;
end


function [miu,sigma] = g1_iter_2d_mu(r,gamma,ModOrder,H_LDPC,rndseed,Nt)
    for k1 = 1:Nt
        rdata = r(k1,:);
        sigmadata = 1./gamma;
        llr_intrlv = qamdemod(rdata.',ModOrder,'OutputType','approxllr','NoiseVariance',sigmadata.');
        llr_intrlv = reshape(llr_intrlv,[],1);
        llr = randdeintrlv(llr_intrlv,rndseed);
        [L0,L1] = converte318(llr);
        [~,llr_ext,~,~] = ldpc_decode_iter(L0,L1,H_LDPC,25);
        llr_ext = randintrlv(real(min(log((1-llr_ext)./(1+llr_ext)+1e-6),15)),rndseed);
        [xhat,vhat] = llr2mv(llr_ext,ModOrder);
        miu(k1,:) = reshape(xhat,size(rdata));
        sigma(k1,:) = reshape(vhat,size(rdata));
    end
    sigma = mean(sigma,1);
end

function [xhat,vhat] = llr2mv(llrin,ModOrder)
    Nbit = log2(ModOrder);
    [c0in,c1in] = converte318(llrin);
    c0_res = reshape(c0in,Nbit,[]);
    c1_res = reshape(c1in,Nbit,[]);
    Bits = (double(dec2bin(0:(ModOrder-1),Nbit))-48).';
    Symbols = qammod(Bits,ModOrder,'InputType','bit','UnitAveragePower',true);
    for k1 = 1:size(c0_res,2)
       Prob_b = repmat(c0_res(:,k1),1,ModOrder);
       t1 = repmat(c1_res(:,k1),1,ModOrder);
       Prob_b(Bits==1) = t1(Bits==1);
       Prob = prod(Prob_b,1);
       p = Prob;
       p = p/sum(p);
       xhat(k1,1) = sum(p.*Symbols);
       vhat(k1,1) = sum(p.*abs(Symbols).^2)-abs(xhat(k1))^2;
    end
end