%%SVD 与GMD基带的方案对比图
%模拟统一SOMP
close all;
clear all;
clc;
%% 与main_irs_upa只有载波数Nc和仿真次数Nsym不同不同
%% system parameter
add_path();
Nsym =20; %number of symbols
%Nsym=100;
N_T =[8 8];N_R=[4 4]; % number of transmit/receive antennas
if length(N_T)==2
    Nt=N_T(1,1)*N_T(1,2);
    Nr=N_R(1,1)*N_R(1,2);
else
    Nt=N_T;
    Nr=N_R;
end
Num_IRS=[16 16];
NumIRS=Num_IRS(1,1)*Num_IRS(1,2);
N_RF = 3; %number of RF chains
N_s = 3;
b = 4; %16QAM,b:modulation order
fc=28e9; % Frequencey
lamada=3e8/fc; % wavelegenth;
SNR_dB=[-10:5:20]; % SNR range
%  SNR_dB=20;
%% channel parameter
los=1;%denote los channel/nlos channel
Nc=64; %subcarrier number
N_cl=8;%number of cluster
N_ray=10;%number of rays per cluster
channel_power = 1;
sigma_2_alpha=channel_power;%cluster power
sigma_ang=7.5;%angle spread
%%
BER_digital_SVD=[];
BER_digital_GMD=[];
BER_hybrid_SOMP_GMD=[];
BER_hybrid_SOMP_SVD=[];
N_tbits=b*Nsym*Nr*Nc; %Total bits
for i_SNR=1:length(SNR_dB)
    SNRdB = SNR_dB(i_SNR);
    N_ebits_digital_SVD = 0;
    N_ebits_digital_GMD = 0;
   N_ebits_hybrid_SOMP_GMD = 0;
   N_ebits_hybrid_SOMP_SVD=0;
    sigma2 = 10^(-SNRdB/10); %generate channel noise convarience
    sigma = sqrt(sigma2);
     tic
    for i_iter=1:Nsym
       
        %% transmit part
        bitseq=randi([0 1],b,N_s);
        symbol_data=bitseq(:)';
        symbol = QAM16_mod(symbol_data,N_s);
        x=symbol.';
        %% mmWave channel matrix generation
        [H1,A_t,A_irs_r] = freqency_sparse_SV_channel0(Nc,N_cl,N_ray, sigma_2_alpha, sigma_ang,Nt,NumIRS,los);
        [H2,A_irs_t,A_r] = freqency_sparse_SV_channel0(Nc,N_cl,N_ray, sigma_2_alpha, sigma_ang,NumIRS,Nr,los);
        %% IRS design
        H_IRS=phase_rotation_design(A_irs_r,A_irs_t,Nc);
        %equivalent channel
        H=zeros(Nr,Nt,Nc);
        for eqv=1:Nc
            H(:,:,eqv)=H2(:,:,eqv)*H_IRS(:,:,eqv)*H1(:,:,eqv);
        end
        %%   RF precoder/combiner
        %% RF PCA
%         [F_RF_PCA,PCAang_F_BB, W_RF_PCA,PCAang_W_BB] = PCAang_pre_and_com(H,N_s,N_RF,N_RF,sigma2,1,0);
        %% RF SOMP
        [F_RF_SOMP,F_BB,W_RF_SOMP,W_BB] = SOMP_pre_and_com(H,A_t, A_r,N_s,N_RF,N_RF,sigma2,1,0);      
        %%  baseband precoding
        for carrier=1:Nc
            [U,S,V] = svd(H(:,:,carrier));
            V_1 = V(:,1:N_s);
            total_power = N_s/sigma2;
            power_allo_equal = eye(N_s)*(total_power/N_s); % equal power allocation
             power_allo_water = water_filling(S,total_power,N_s);
           %%  digital SVD
           precode_symbol_digital_SVD = V_1*power_allo_water*x*(N_s/total_power);
           %% hybrid SVD
           F_hybrid_SVD_SOMP=sqrt(N_s/N_RF)*F_RF_SOMP*F_BB(:,:,carrier);
           precode_symbol_hybrid_SVD_SOMP = F_hybrid_SVD_SOMP*power_allo_water*x*(N_s/total_power);
            %% digital GMD
            [G,M,D] = gmd(U,S(:,1:N_s),V(:,1:N_s));
            precode_symbol_digital_GMD = D*power_allo_equal*x*(N_s/total_power);%% full digital
            %% hybrid GMD
            [U2,S2,V2] = svd(W_RF_SOMP'*H(:,:,carrier)*F_RF_SOMP);
            [G_SOMP,M2,D2] = gmd(U2,S2,V2);
            F_hybrid_GMD_SOMP = sqrt(N_s/N_RF)*F_RF_SOMP*D2;
            precode_symbol_hybrid_GMD_SOMP = F_hybrid_GMD_SOMP*power_allo_equal*x*(N_s/total_power);

            %% receive part
            noise = sigma*(randn(Nr,1)+1i*randn(Nr,1));  % generate the noise
            y_digital_SVD=H(:,:,carrier)* precode_symbol_digital_SVD+noise;
            y_digital_GMD = H(:,:,carrier)*precode_symbol_digital_GMD + noise;
            y_hybrid_SOMP_SVD = H(:,:,carrier)*  precode_symbol_hybrid_SVD_SOMP + noise;
            y_hybrid_SOMP_GMD = H(:,:,carrier)*  precode_symbol_hybrid_GMD_SOMP + noise;
             

            %% combining
             GH = G';UH=U';
            G_1 = GH(1:N_s,:);
            
            y_hat_digital_SVD=UH(1:N_s,:)*y_digital_SVD;
            
            y_hat_digital_GMD = G_1*y_digital_GMD;
                    
            W_hybrid_SVD_SOMP= W_BB(:,:,carrier)'*W_RF_SOMP';
            y_hat_hybrid_SOMP_SVD = sqrt(N_s/N_RF)*W_hybrid_SVD_SOMP*y_hybrid_SOMP_SVD;
            
            W_hybrid_GMD_SOMP= G_SOMP'*W_RF_SOMP';
            y_hat_hybrid_SOMP_GMD = sqrt(N_s/N_RF)*W_hybrid_GMD_SOMP*y_hybrid_SOMP_GMD;
            
            
            %%%%%%%%%%%%%%%%%%%%% SVD Decoder %%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
              symbol_sliced_digital_SVD = VBLAST_decoder(y_hat_digital_SVD,N_s,S(:,1:N_s));
              demapper_data_digital_SVD = QAM16_demapper(symbol_sliced_digital_SVD);
            
             symbol_sliced_hybrid_SOMP_SVD = VBLAST_decoder(y_hat_hybrid_SOMP_SVD,N_s,S(:,1:N_s));
             demapper_data_hybrid_SOMP_SVD = QAM16_demapper(symbol_sliced_hybrid_SOMP_SVD);

            %%%%%%%%%%%%%%%%%%%% GMD-SIC decoder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            symbol_sliced_digital_GMD = VBLAST_decoder(y_hat_digital_GMD,N_s,M(:,1:N_s));
            demapper_data_digital_GMD = QAM16_demapper(symbol_sliced_digital_GMD);

            symbol_sliced_hybrid_SOMP_GMD = VBLAST_decoder(y_hat_hybrid_SOMP_GMD,N_s,M2);
            demapper_data_hybrid_SOMP_GMD = QAM16_demapper(symbol_sliced_hybrid_SOMP_GMD);
            

            %%%%%%%%%%%%%%%%%%% error counting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            N_ebits_digital_SVD = N_ebits_digital_SVD + sum(symbol_data~= demapper_data_digital_SVD);         
            N_ebits_hybrid_SOMP_SVD = N_ebits_hybrid_SOMP_SVD + sum(symbol_data~= demapper_data_hybrid_SOMP_SVD);                
            N_ebits_digital_GMD = N_ebits_digital_GMD + sum(symbol_data~= demapper_data_digital_GMD);
            N_ebits_hybrid_SOMP_GMD = N_ebits_hybrid_SOMP_GMD + sum(symbol_data~= demapper_data_hybrid_SOMP_GMD);
           
        end
       
       fprintf('SNR=%d dB,iter=%d\n',SNR_dB(i_SNR),i_iter );
    end
      toc
    BER_digital_SVD(i_SNR) = N_ebits_digital_SVD/N_tbits;
    BER_digital_GMD(i_SNR) = N_ebits_digital_GMD/N_tbits;
    BER_hybrid_SOMP_GMD(i_SNR) = N_ebits_hybrid_SOMP_GMD/N_tbits;
    BER_hybrid_SOMP_SVD(i_SNR) = N_ebits_hybrid_SOMP_SVD/N_tbits;

end
%%  BER PLOT

semilogy(SNR_dB,smooth(BER_digital_SVD),'b-.','Linewidth',1.5)
hold on
semilogy(SNR_dB,smooth( BER_hybrid_SOMP_SVD),'r-s','Linewidth',1.5)
hold on
semilogy(SNR_dB,smooth(BER_digital_GMD),'k-.','Linewidth',1.5)
hold on
semilogy(SNR_dB,smooth( BER_hybrid_SOMP_GMD),'g-s','Linewidth',1.5)
hold on

xlabel('SNR (dB)')
ylabel('BER')
legend('Fully digital SVD','Hybrid-SVD','Fully digital GMD','Hybrid-GMD');
grid on
ylim([10^-5 10^-1])

