close all;
clear all;
clc;
%% system parameter
add_path();
Nsym =200; %number of symbols
N_T =[8 8];N_R=[4 4]; % number of transmit/receive antennas
Num_IRS=[16 16];
NumIRS=Num_IRS(1,1)*Num_IRS(1,2);

if length(N_T)==2
    Nt=N_T(1,1)*N_T(1,2);
    Nr=N_R(1,1)*N_R(1,2);
else
    Nt=N_T;
    Nr=N_R;
end
N_RF = 3; %number of RF chains
N_s = 3;
b = 4; %16QAM,b:modulation order
fc=28e9; % Frequencey
lamada=3e8/fc; % wavelegenth;
SNR_dB=[-10:5:20]; % SNR range
%  SNR_dB=25;
%% channel parameter
los=1;%denote los channel/nlos channel
Nc=64; %subcarrier number
N_cl=8;%number of cluster
N_ray=10;%number of rays per cluster
channel_power = 1;
sigma_2_alpha=channel_power;%cluster power
sigma_ang=7.5;%angle spread
%% codebook generate
rho1=1;%over samping
rho2=2;
rho3=3;
Q_A=[log2(N_T(1)) log2(N_R(1))];
DFT_set_r1 = Q_codebook(N_R, Q_A(2), rho1); DFT_set_t1 = Q_codebook(N_T,Q_A(1), rho1);
DFT_set_r2 = Q_codebook(N_R,  Q_A(2) ,rho2); DFT_set_t2 = Q_codebook(N_T, Q_A(1),rho2);
DFT_set_r3 = Q_codebook(N_R,  Q_A(2), rho3); DFT_set_t3 = Q_codebook(N_T, Q_A(1),rho3);
%%
BER_hybrid_PCA = [];BER_hybrid_SOMP = [];BER_hybrid_MAO=[];BER_digital_GMD=[];
BER_hybrid_SOMP1 = []; BER_hybrid_SOMP2 = [];BER_hybrid_SOMP3 = [];
N_tbits=b*Nsym*Nr*Nc; %Total bits
for i_SNR=1:length(SNR_dB)
    SNRdB = SNR_dB(i_SNR);
    N_ebits_hybrid_PCA = 0;N_ebits_hybrid_SOMP = 0;N_ebits_hybrid_MAO = 0;N_ebits_digital_GMD = 0;
    N_ebits_hybrid_SOMP1 = 0;N_ebits_hybrid_SOMP2 = 0;N_ebits_hybrid_SOMP3 = 0;
    sigma2 = 10^(-SNRdB/10); %generate channel noise convarience
    sigma = sqrt(sigma2);  
    
    for i_iter=1:Nsym
         tic
        %% transmit part
        bitseq=randi([0 1],b,N_s);
        symbol_data=bitseq(:)';
        symbol = QAM16_mod(symbol_data,N_s);
        x=symbol.';
        %% mmWave channel matrix generation
%         [H, A_BS,A_MS] = channel_f(Nt, Nr, 1, Nc);
%       这里的信道根据Los/Nlos径略有修改信道模型的line86-92
        [H1,A_t,A_irs_r] = freqency_sparse_SV_channel0(Nc,N_cl,N_ray, sigma_2_alpha, sigma_ang,N_T,Num_IRS,los);
        [H2,A_irs_t,A_r] = freqency_sparse_SV_channel0(Nc,N_cl,N_ray, sigma_2_alpha, sigma_ang,Num_IRS,N_R,los);
        %% IRS design
        H_IRS=phase_rotation_design(A_irs_r,A_irs_t,Nc);
        %equivalent channel
        H=zeros(Nr,Nt,Nc);
        for eqv=1:Nc
            H(:,:,eqv)=H2(:,:,eqv)*H_IRS(:,:,eqv)*H1(:,:,eqv);
%           H(:,:,eqv)=H2(:,:,eqv)*H1(:,:,eqv);
        end
        %%   RF precoder/combiner
        %% RF PCA
       [F_RF_PCA,PCAang_F_BB, W_RF_PCA,PCAang_W_BB] = PCAang_pre_and_com(H,N_s,N_RF,N_RF,sigma2,1,0);
        %% RF SOMP
        [F_RF_SOMP,~,W_RF_SOMP,~] = SOMP_pre_and_com(H,A_t, A_r,N_s,N_RF,N_RF,sigma2,1,0);   
        [F_RF_SOMP_1,~,W_RF_SOMP_1,~] = SOMP_pre_and_com(H,DFT_set_t1, DFT_set_r1,N_s,N_RF,N_RF,sigma2,1,0);
        [F_RF_SOMP_2,~,W_RF_SOMP_2,~] = SOMP_pre_and_com(H,DFT_set_t2, DFT_set_r2,N_s,N_RF,N_RF,sigma2,1,0);
        [F_RF_SOMP_3,~,W_RF_SOMP_3,~] = SOMP_pre_and_com(H,DFT_set_t3, DFT_set_r3,N_s,N_RF,N_RF,sigma2,1,0);
         

        %%  baseband precoding
        for carrier=1:Nc
            [U,S,V] = svd(H(:,:,carrier));
            V_1 = V(:,1:N_s);
            total_power = N_s/sigma2;
            power_allo_equal = eye(N_s)*(total_power/N_s); % equal power allocation
            %% digital GMD
            [G,M,D] = gmd(U,S(:,1:N_s),V(:,1:N_s));
            precode_symbol_digital_GMD = D*power_allo_equal*x*(N_s/total_power);%% full digital
            %% hybrid GMD
            [U1,S1,V1] = svd(W_RF_PCA'*H(:,:,carrier)*F_RF_PCA);
            [G_PCA,M1,D1] = gmd(U1,S1,V1);
            F_hybrid_GMD_PCA = sqrt(N_s/N_RF)*F_RF_PCA*D1;
            precode_symbol_hybrid_GMD_PCA = F_hybrid_GMD_PCA*power_allo_equal*x*(N_s/total_power);
            
            [U2,S2,V2] = svd(W_RF_SOMP'*H(:,:,carrier)*F_RF_SOMP);
            [G_SOMP,M2,D2] = gmd(U2,S2,V2);
            F_hybrid_GMD_SOMP = sqrt(N_s/N_RF)*F_RF_SOMP*D2;
            precode_symbol_hybrid_GMD_SOMP = F_hybrid_GMD_SOMP*power_allo_equal*x*(N_s/total_power);
            
            [U2,S2,V2] = svd(W_RF_SOMP_1'*H(:,:,carrier)*F_RF_SOMP_1);
            [G_SOMP1,M21,D2] = gmd(U2,S2,V2);
            F_hybrid_GMD_SOMP1 = sqrt(N_s/N_RF)*F_RF_SOMP_1*D2;
            precode_symbol_hybrid_GMD_SOMP1 = F_hybrid_GMD_SOMP1*power_allo_equal*x*(N_s/total_power);
            
             [U2,S2,V2] = svd(W_RF_SOMP_2'*H(:,:,carrier)*F_RF_SOMP_2);
            [G_SOMP2,M22,D2] = gmd(U2,S2,V2);
            F_hybrid_GMD_SOMP2 = sqrt(N_s/N_RF)*F_RF_SOMP_2*D2;
            precode_symbol_hybrid_GMD_SOMP2 = F_hybrid_GMD_SOMP2*power_allo_equal*x*(N_s/total_power);
            
            [U2,S2,V2] = svd(W_RF_SOMP_3'*H(:,:,carrier)*F_RF_SOMP_3);
            [G_SOMP3,M23,D2] = gmd(U2,S2,V2);
            F_hybrid_GMD_SOMP3 = sqrt(N_s/N_RF)*F_RF_SOMP_3*D2;
            precode_symbol_hybrid_GMD_SOMP3 = F_hybrid_GMD_SOMP3*power_allo_equal*x*(N_s/total_power);
            
            
%             
%             [U3,S3,V3] = svd(W_RF_MAO'*H(:,:,carrier)*F_RF_MAO);
%             [G_MAO,M3,D3] = gmd(U3,S3,V3);
%             F_hybrid_GMD_MAO = sqrt(N_s/N_RF)*F_RF_MAO*D3;
%             precode_symbol_hybrid_GMD_MAO = F_hybrid_GMD_MAO*power_allo_equal*x*(N_s/total_power);
            %% receive part
            noise = sigma*(randn(Nr,1)+1i*randn(Nr,1));  % generate the noise
            y_digital_GMD = H(:,:,carrier)*precode_symbol_digital_GMD + noise;
            y_hybrid_PCA = H(:,:,carrier)*  precode_symbol_hybrid_GMD_PCA + noise;
            y_hybrid_SOMP = H(:,:,carrier)*  precode_symbol_hybrid_GMD_SOMP + noise;
            y_hybrid_SOMP1 = H(:,:,carrier)*  precode_symbol_hybrid_GMD_SOMP1 + noise;
            y_hybrid_SOMP2 = H(:,:,carrier)*  precode_symbol_hybrid_GMD_SOMP2 + noise;
            y_hybrid_SOMP3 = H(:,:,carrier)*  precode_symbol_hybrid_GMD_SOMP3 + noise;
%             y_hybrid_MAO = H(:,:,carrier)*  precode_symbol_hybrid_GMD_MAO + noise;
            
            %% combining
            
            % digital GMD
            GH = G';
            G_1 = GH(1:N_s,:);
            y_hat_digital_GMD = G_1*y_digital_GMD;
            
            % hybrid GMD
            W_hybrid_GMD_PCA = G_PCA'*W_RF_PCA';
            y_hat_hybrid_PCA = sqrt(N_s/N_RF)*W_hybrid_GMD_PCA*y_hybrid_PCA;
            
            W_hybrid_GMD_SOMP= G_SOMP'*W_RF_SOMP';
            y_hat_hybrid_SOMP = sqrt(N_s/N_RF)*W_hybrid_GMD_SOMP*y_hybrid_SOMP;
            
            W_hybrid_GMD_SOMP1= G_SOMP1'*W_RF_SOMP_1';
            y_hat_hybrid_SOMP1 = sqrt(N_s/N_RF)*W_hybrid_GMD_SOMP1*y_hybrid_SOMP1;
            
            W_hybrid_GMD_SOMP2= G_SOMP2'*W_RF_SOMP_2';
            y_hat_hybrid_SOMP2 = sqrt(N_s/N_RF)*W_hybrid_GMD_SOMP2*y_hybrid_SOMP2;
            
            W_hybrid_GMD_SOMP3= G_SOMP3'*W_RF_SOMP_3';
            y_hat_hybrid_SOMP3 = sqrt(N_s/N_RF)*W_hybrid_GMD_SOMP3*y_hybrid_SOMP3;
%             
%             W_hybrid_GMD_MAO = G_MAO'*W_RF_MAO';
%             y_hat_hybrid_MAO = sqrt(N_s/N_RF)*W_hybrid_GMD_MAO*y_hybrid_MAO;
            %%%%%%%%%%%%%%%%%%%% GMD-SIC decoder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            symbol_sliced_digital_GMD = VBLAST_decoder(y_hat_digital_GMD,N_s,M(:,1:N_s));
            demapper_data_digital_GMD = QAM16_demapper(symbol_sliced_digital_GMD);
            symbol_sliced_hybrid_PCA = VBLAST_decoder(y_hat_hybrid_PCA,N_s,M1);
            demapper_data_hybrid_PCA = QAM16_demapper(symbol_sliced_hybrid_PCA);
            symbol_sliced_hybrid_SOMP = VBLAST_decoder(y_hat_hybrid_SOMP,N_s,M2);
            demapper_data_hybrid_SOMP = QAM16_demapper(symbol_sliced_hybrid_SOMP);
            
              symbol_sliced_hybrid_SOMP1 = VBLAST_decoder(y_hat_hybrid_SOMP1,N_s,M21);
            demapper_data_hybrid_SOMP1 = QAM16_demapper(symbol_sliced_hybrid_SOMP1);
              symbol_sliced_hybrid_SOMP2 = VBLAST_decoder(y_hat_hybrid_SOMP2,N_s,M22);
            demapper_data_hybrid_SOMP2 = QAM16_demapper(symbol_sliced_hybrid_SOMP2);
              symbol_sliced_hybrid_SOMP3 = VBLAST_decoder(y_hat_hybrid_SOMP3,N_s,M23);
            demapper_data_hybrid_SOMP3 = QAM16_demapper(symbol_sliced_hybrid_SOMP3);
                
%             symbol_sliced_hybrid_MAO = VBLAST_decoder(y_hat_hybrid_MAO,N_s,M3);
%             demapper_data_hybrid_MAO = QAM16_demapper(symbol_sliced_hybrid_MAO);
            %%%%%%%%%%%%%%%%%%% error counting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            N_ebits_digital_GMD = N_ebits_digital_GMD + sum(symbol_data~= demapper_data_digital_GMD);
            N_ebits_hybrid_PCA = N_ebits_hybrid_PCA + sum(symbol_data~= demapper_data_hybrid_PCA);
            N_ebits_hybrid_SOMP = N_ebits_hybrid_SOMP + sum(symbol_data~= demapper_data_hybrid_SOMP);
            N_ebits_hybrid_SOMP1 = N_ebits_hybrid_SOMP1 + sum(symbol_data~= demapper_data_hybrid_SOMP1);
            N_ebits_hybrid_SOMP2 = N_ebits_hybrid_SOMP2 + sum(symbol_data~= demapper_data_hybrid_SOMP2);
            N_ebits_hybrid_SOMP3 = N_ebits_hybrid_SOMP3 + sum(symbol_data~= demapper_data_hybrid_SOMP3);
%             N_ebits_hybrid_MAO = N_ebits_hybrid_MAO + sum(symbol_data~= demapper_data_hybrid_MAO);
        end
       fprintf('SNR=%d ,iter=%d\n',SNR_dB(i_SNR),i_iter );
     toc 
    end
    BER_digital_GMD(i_SNR) = N_ebits_digital_GMD/N_tbits;
    BER_hybrid_PCA(i_SNR) = N_ebits_hybrid_PCA/N_tbits;
    BER_hybrid_SOMP(i_SNR) = N_ebits_hybrid_SOMP/N_tbits;
    BER_hybrid_SOMP1(i_SNR) = N_ebits_hybrid_SOMP1/N_tbits;
    BER_hybrid_SOMP2(i_SNR) = N_ebits_hybrid_SOMP2/N_tbits;
    BER_hybrid_SOMP3(i_SNR) = N_ebits_hybrid_SOMP3/N_tbits;
%     BER_hybrid_MAO(i_SNR) = N_ebits_hybrid_MAO/N_tbits;
end
semilogy(SNR_dB,smooth(BER_digital_GMD),'k-.','Linewidth',1.5)
hold on
semilogy(SNR_dB,smooth( BER_hybrid_PCA),'b-o','Linewidth',1.5)
hold on

semilogy(SNR_dB,smooth( BER_hybrid_SOMP),'g-s','Linewidth',1.5)
hold on
semilogy(SNR_dB,smooth( BER_hybrid_SOMP1),'r-v','Linewidth',1.5)
hold on
semilogy(SNR_dB,smooth( BER_hybrid_SOMP2),'r-+','Linewidth',1.5)
hold on
semilogy(SNR_dB,smooth( BER_hybrid_SOMP3),'r-d','Linewidth',1.5)
hold on
% 
% semilogy(SNR_dB,smooth( BER_hybrid_MAO),'c-*','Linewidth',1.5)

xlabel('SNR (dB)')
ylabel('BER')
legend('Full-digital','PCA','Ideal-SOMP','\rho=1 SOMP','\rho=2 SOMP','\rho=3 SOMP');
grid on
ylim([10^-5 10^-1])

