clear,clc;
Ka = 100;
Nr = 100;Nspath = 18;
PilotBits = 7;PilotNum = 2^PilotBits;
P = 5;N = 3;
L = 128;
ModOrder = 2;ModBits = log2(ModOrder);
Lall = 640;Ldata = Lall-L;
Nbits = 100;
DataBits = Nbits-PilotBits;
% CRCBits = 11;CRCstr = num2str(CRCBits);
% PCBits = Ldata*ModBits-DataBits-CRCBits;

PCBits = Ldata*ModBits-DataBits;
s = rng;
H_LDPC = ldpc_generate(PCBits,Ldata*ModBits,3,2,347234);
[H_LDPC, G_LDPC] = ldpc_h2g(H_LDPC);
rndseed = 774;
rng(s);


indpilot = ceil(linspace(1,Lall,L));
inddata = 1:Lall;inddata(indpilot) = [];
% Transform Matrix construction
Fmat = dftmtx(Lall+P);
Fmat_pilot = Fmat(P+indpilot,1:P);
Fmat_data = Fmat(P+inddata,1:P);
Fmat_spatial = dftmtx(Nr)/sqrt(Nspath);% Ensure normalized power of channel mat.

s = RandStream('mt19937ar','Seed',15,'NormalTransform','Polar');
A = 1/sqrt(2)*(randn(s,L,PilotNum)+1j*randn(s,L,PilotNum));
A_mat_f = [];

for k1 = 1:PilotNum
    A_mat_f = [A_mat_f, sparse(diag(A(:,k1)))*Fmat_pilot];
end


SNR = (-9:0.3:-6);

N0 = 1./(10.^(SNR/10));
Eb_N0 = 10*log10(Lall./(Nbits*N0))
SimNum = 60;

for ksnr = 1:length(N0)
    tic;
    for ksim = 1:SimNum
    % User data generation and transmission

        Hact = randperm(PilotNum,Ka);
        Hact_st = zeros(PilotNum,1);
        H_sparse_stitched = zeros(PilotNum*P,Nr);
        Z = zeros(Lall,Nr);
        Data_raw_p1 = zeros(DataBits,Ka);
        Data_raw_p2 = zeros(1,Ka);
        Data_Coded = zeros(Ldata*ModBits,Ka);
        ModSym = zeros(Ldata,Ka);
        PilotSym = zeros(L,Ka);
        Data = zeros(Lall,Ka);
        H_sparse = cell(1,Ka);
        H_all = cell(1,Ka);
        H_mask_delay = zeros(P,Nr);
        for k1 = 1:Ka
            Hact_st(Hact(k1)) = 1;
            Data_raw_p1(:,k1) = randi([0,1],DataBits,1);
            % Data_raw_p1_CRC(:,k1) = a5gCRCEncode(Data_raw_p1(:,k1).', CRCstr).';
            Data_raw_p2(k1) = randi([1,PilotNum]);
            while(ismember(Data_raw_p2(k1),Data_raw_p2(1:k1-1)))
                Data_raw_p2(k1) = randi([1,PilotNum]);
            end
            Data_Coded(:,k1) = randintrlv(mod(Data_raw_p1(:,k1).'*G_LDPC,2).',rndseed);
            % Data_Coded_temp = reshape(Data_Coded(:,k1),ModBits,[]);
            ModSym (:,k1) = qammod(Data_Coded(:,k1),ModOrder,"InputType","bit","UnitAveragePower",true);
            PilotSym (:,k1) = A(:,Data_raw_p2(k1));
            Data(indpilot,k1) = PilotSym (:,k1);
            Data(inddata,k1) = ModSym (:,k1);


            H_temp = 1/sqrt(2)*(randn(P,Nr)+1j*randn(P,Nr));

            for k2 = 1:Nr
                H_mask_delay(:,k2) = randperm(P).'>(P-N);
            end
            H_mask_ant = repmat(randperm(Nr)>(Nr-Nspath),P,1);
            H_mask = H_mask_delay.*H_mask_ant;
            H_sparse{k1} = H_mask.*H_temp/sqrt(N);
            H_sparse_stitched((Data_raw_p2(k1)-1)*P+1:(Data_raw_p2(k1))*P,:) = H_sparse{k1}; 
            H_pilot = Fmat_pilot*H_sparse{k1}*Fmat_spatial';
            H_data = Fmat_data*H_sparse{k1}*Fmat_spatial';
            H_all{k1}(indpilot,:) = H_pilot;
            H_all{k1}(inddata,:) = H_data;

            Z = Z+H_all{k1}.*repmat(Data(:,k1),1,Nr);

        end
        Noise = sqrt(N0(ksnr)/2)*(randn(Lall,Nr)+1j*randn(Lall,Nr));
        Y = Z+Noise;
        Y_SIC = Y;
        % SIC-VAMP detection
        SICiter = 3;

        maflg = ones(1,Ka);
        for ksic = 1:SICiter
            Y_pilot = Y_SIC(indpilot,:);
            Y_pilot_angle = Y_pilot*Fmat_spatial*Nspath/Nr;

            % Sparse domain Channel Estimation
            % [H_est,tau] = EM_VAMP_MMV_2d_SVD_c(Y_pilot_angle,A_mat_f,N0,N/P*Ka/PilotNum,Nspath/Nr,1,H_all);
            [H_est,tau] = VAMP_MMV_2d_SVD_c(Y_pilot_angle,A_mat_f,N0(ksnr),N/P*Ka/PilotNum,Nspath/Nr,1/sqrt(N));
            % stem3(abs(H_est),'b*');hold on;
            % stem3(abs(H_sparse_stitched),'ro');
            % legend('VAMP','original')
            MSE_VAMP = mean(abs(H_est(H_sparse_stitched>0)-H_sparse_stitched(H_sparse_stitched>0)).^2)
            % MSE_VAMP = mean(abs(H_est-H_sparse_stitched).^2,'all')
            kuser = 1;
            Pilotidx = zeros(1,PilotNum);
            H_est_all = cell(1,PilotNum);
            H_est_data_freq_spatial = cell(1,PilotNum);
            for k1 = 1:PilotNum
                if max(abs(H_est((k1-1)*P+1:k1*P,:)),[],'all')>0.01
                    H_est_all{kuser} = H_est((k1-1)*P+1:k1*P,:);
                    H_est_data_freq_spatial{kuser} = Fmat_data*H_est((k1-1)*P+1:k1*P,:)*Fmat_spatial';
                    Pilotidx(kuser) = k1;
                    kuser = kuser + 1;
                else
                    continue;
                end
            end
            Nuser_est = kuser-1;
            if Nuser_est == 0
                break;
            end
            H_reorg = zeros(Nr,Nuser_est,Ldata);
            for k1 = 1:Ldata
                for k2 = 1:Nuser_est
                    H_reorg(:,k2,k1) = H_est_data_freq_spatial{k2}(k1,:).';
                end
            end

            % Detection
            Y_data = Y_SIC(inddata,:);
            [x,~] = VAMP_2d_svd_iter_c(Y_data.',H_reorg,N0(ksnr),ModOrder,H_LDPC,rndseed);
            Bitdec = zeros(DataBits,Nuser_est);
            for k1 = 1:Nuser_est
                llr_intrlv = qamdemod(x(k1,:),ModOrder,'OutputType','approxllr','NoiseVariance',N0(ksnr));
                llr_intrlv = reshape(llr_intrlv,[],1);
                llr = randdeintrlv(llr_intrlv,rndseed);
                [L0,L1] = converte318(llr);
                [llr1_out,~,correctflg,~] = ldpc_decode_iter(L0,L1,H_LDPC,25);
                Bitdec(:,k1) = llr1_out>0;
                % Bitchk = a5gCRCEncode(Bitdec(1:DataBits,k1).',CRCstr).';
                if correctflg
                    H_rec = zeros(Lall,Nr);
                    H_data_rec = Fmat_data*H_est_all{k1}*Fmat_spatial';
                    H_pilot_rec = Fmat_pilot*H_est_all{k1}*Fmat_spatial';
                    H_rec(indpilot,:) = H_pilot_rec;
                    H_rec(inddata,:) = H_data_rec;
                    Bit_rec_Coded = randintrlv(mod(Bitdec(:,k1).'*G_LDPC,2).',rndseed);
                    ModSym_rec = qammod(Bit_rec_Coded,ModOrder,'InputType','bit',"UnitAveragePower",true);
                    PilotSym_rec = A(:,Pilotidx(k1));

                    Data_rec = zeros(Lall,1);

                    Data_rec(indpilot,1) = PilotSym_rec;
                    Data_rec(inddata,1) = ModSym_rec;
                    Z_rec = H_rec.*repmat(Data_rec,1,Nr);
                    Y_SIC = Y_SIC-Z_rec;
                end
            end

            for k1 = 1:Nuser_est
                for k2 = 1:Ka
                    if sum(mod(Bitdec(1:DataBits,k1)+Data_raw_p1(:,k2),2))==0
                        maflg(k2) = 0;
                        break;
                    end
                end
            end

            if sum(maflg) == 0
                break;
            end
        end
        manum(ksim) = sum(maflg);
    end
    toc;

    PUPE(ksnr) = sum(manum)/(SimNum*Ka)
end
