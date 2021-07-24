
close all;
clear all;
clc;
%% 与main_irs_upa只有载波数Nc和仿真次数Nsym不同不同
%% system parameter
add_path();
Nsym =300; %number of symbols
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
SNR_dB=[-10:5:30]; % SNR range
%  SNR_dB=20;
%% channel parameter
los=1;%denote los channel/nlos channel
Nc=3; %subcarrier number
N_cl=8;%number of cluster
N_ray=10;%number of rays per cluster
channel_power = 1;
sigma_2_alpha=channel_power;%cluster power
sigma_ang=7.5;%angle spread
%% codebook generate
rho1=1;%over samping
rho2=2;
rho3=2;
Q_A=[log2(N_T(1)) log2(N_R(1))];
DFT_set_r1 = Q_codebook(N_R, Q_A(2), rho1); DFT_set_t1 = Q_codebook(N_T,Q_A(1), rho1);
DFT_set_r2 = Q_codebook(N_R,  Q_A(2) ,rho2); DFT_set_t2 = Q_codebook(N_T, Q_A(1),rho2);
DFT_set_r3 = Q_codebook(N_R,  Q_A(2), rho3); DFT_set_t3 = Q_codebook(N_T, Q_A(1),rho3);
%%
N_tbits=b*Nsym*Nr*Nc; %Total bits


R_hybrid_PCA=zeros(1,length(SNR_dB));
R_hybrid_SOMP=zeros(1,length(SNR_dB));
R_hybrid_SOMP1=zeros(1,length(SNR_dB));
R_hybrid_SOMP2=zeros(1,length(SNR_dB));
R_hybrid_SOMP3=zeros(1,length(SNR_dB));
R_hybrid_SOMP4=zeros(1,length(SNR_dB));
R_hybrid_SOMP5=zeros(1,length(SNR_dB));
R_hybrid_SOMP6=zeros(1,length(SNR_dB));
R_hybrid_SOMP7=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB)
    SNRdB = SNR_dB(i);
    snr =10.^(SNR_dB(i) / 10);
 
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
%         [H, A_BS,A_MS] = channel_f(Nt, Nr, 1, Nc);
%       这里的信道根据Los/Nlos径略有修改信道模型的line86-92
        [H1,A_t,A_irs_r] = freqency_sparse_SV_channel0(Nc,N_cl,N_ray, sigma_2_alpha, sigma_ang,Nt,NumIRS,los);
        [H2,A_irs_t,A_r] = freqency_sparse_SV_channel0(Nc,N_cl,N_ray, sigma_2_alpha, sigma_ang,NumIRS,Nr,los);
        [H3,A_tt,A_rr] = freqency_sparse_SV_channel0(Nc,N_cl,N_ray, sigma_2_alpha, sigma_ang,Nt,Nr,0);
        %% IRS design
        H_IRS=phase_rotation_design(A_irs_r,A_irs_t,Nc);
        %equivalent channel
        H=zeros(Nr,Nt,Nc); H_temp=zeros(Nr,Nt,Nc); 
        for eqv=1:Nc

            H(:,:,eqv)=H2(:,:,eqv)*H_IRS(:,:,eqv)*H1(:,:,eqv);
            H_temp(:,:,eqv)=H3(:,:,eqv);
        end
        %%   RF precoder/combiner
        %% RF PCA
        [F_RF_PCA,PCAang_F_BB, W_RF_PCA,PCAang_W_BB] = PCAang_pre_and_com(H,N_s,N_RF,N_RF,sigma2,1,0);
        %% RF SOMP
        %ideal
        [F_RF_SOMP,~,W_RF_SOMP,~] = SOMP_pre_and_com(H,A_t, A_r,N_s,N_RF,N_RF,sigma2,1,0);   
        [F_RF_SOMP_1,~,W_RF_SOMP_1,~] = SOMP_pre_and_com(H_temp,A_tt, A_rr,N_s,N_RF,N_RF,sigma2,1,0);
        %codebook
        %rho=1;
        [F_RF_SOMP_2,~,W_RF_SOMP_2,~] = SOMP_pre_and_com(H,DFT_set_t1, DFT_set_r1,N_s,N_RF,N_RF,sigma2,1,0);
        [F_RF_SOMP_3,~,W_RF_SOMP_3,~] = SOMP_pre_and_com(H_temp,DFT_set_t1, DFT_set_r1,N_s,N_RF,N_RF,sigma2,1,0);
        %rho=2;
        [F_RF_SOMP_4,~,W_RF_SOMP_4,~] = SOMP_pre_and_com(H,DFT_set_t2, DFT_set_r2,N_s,N_RF,N_RF,sigma2,1,0);
        [F_RF_SOMP_5,~,W_RF_SOMP_5,~] = SOMP_pre_and_com(H_temp,DFT_set_t2, DFT_set_r2,N_s,N_RF,N_RF,sigma2,1,0);
        
        %rho=3
        [F_RF_SOMP_6,~,W_RF_SOMP_6,~] = SOMP_pre_and_com(H,DFT_set_t3, DFT_set_r3,N_s,N_RF,N_RF,sigma2,1,0);
        [F_RF_SOMP_7,~,W_RF_SOMP_7,~] = SOMP_pre_and_com(H_temp,DFT_set_t3, DFT_set_r3,N_s,N_RF,N_RF,sigma2,1,0);
        %%  baseband precoding
        for carrier=1:Nc
            [U,S,V] = svd(H(:,:,carrier));
            V_1 = V(:,1:N_s);
            total_power = N_s/sigma2;
            power_allo_equal = eye(N_s)*(total_power/N_s); % equal power allocation
            %% digital GMD
            [G(:,:,carrier),M,D(:,:,carrier)] = gmd(U,S(:,1:N_s),V(:,1:N_s));
             GH(:,:,carrier) = (G(:,:,carrier))';
             QQ=(GH(:,:,carrier));
             G_1(:,:,carrier) = QQ(1:N_s,:);
            %% hybrid GMD
            [U1,S1,V1] = svd(W_RF_PCA'*H(:,:,carrier)*F_RF_PCA);
            [G_PCA(:,:,carrier),M1,D1(:,:,carrier)] = gmd(U1,S1,V1);
         
            
            [U2,S2,V2] = svd(W_RF_SOMP'*H(:,:,carrier)*F_RF_SOMP);
            [G_SOMP(:,:,carrier),M2,D2(:,:,carrier)] = gmd(U2,S2,V2);
       
            
            [U2,S2,V2] = svd(W_RF_SOMP_1'*H_temp(:,:,carrier)*F_RF_SOMP_1);
            [G_SOMP1(:,:,carrier),M21,D21(:,:,carrier)] = gmd(U2,S2,V2);
         
             [U2,S2,V2] = svd(W_RF_SOMP_2'*H(:,:,carrier)*F_RF_SOMP_2);
            [G_SOMP2(:,:,carrier),M22,D22(:,:,carrier)] = gmd(U2,S2,V2);
 
            [U2,S2,V2] = svd(W_RF_SOMP_3'*H_temp(:,:,carrier)*F_RF_SOMP_3);
            [G_SOMP3(:,:,carrier),M23,D23(:,:,carrier)] = gmd(U2,S2,V2);
% 

             [U2,S2,V2] = svd(W_RF_SOMP_4'*H(:,:,carrier)*F_RF_SOMP_4);
            [G_SOMP4(:,:,carrier),M24,D24(:,:,carrier)] = gmd(U2,S2,V2);
 
            [U2,S2,V2] = svd(W_RF_SOMP_5'*H_temp(:,:,carrier)*F_RF_SOMP_5);
            [G_SOMP5(:,:,carrier),M25,D25(:,:,carrier)] = gmd(U2,S2,V2);
            
            [U2,S2,V2] = svd(W_RF_SOMP_6'*H(:,:,carrier)*F_RF_SOMP_6);
            [G_SOMP6(:,:,carrier),M26,D26(:,:,carrier)] = gmd(U2,S2,V2);
 
            [U2,S2,V2] = svd(W_RF_SOMP_7'*H_temp(:,:,carrier)*F_RF_SOMP_7);
            [G_SOMP7(:,:,carrier),M27,D27(:,:,carrier)] = gmd(U2,S2,V2);
% 
% 
%             [U3,S3,V3] = svd(W_RF_MAO'*H(:,:,carrier)*F_RF_MAO);
%             [G_MAO(:,:,carrier),M3,D3(:,:,carrier)] = gmd(U3,S3,V3);
        end       
          %%%%%%%%%%%%%%%%%%%%spectrum efficiency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             R_digital_GMD(:,i) =  R_digital_GMD(:,i)+SE_BER(1, H, G_1,0, D, 0,channel_power, snr);
            R_hybrid_PCA(:,i)= R_hybrid_PCA(:,i)+SE_BER(2, H, G_PCA, W_RF_PCA, D1, F_RF_PCA, channel_power, snr);
            R_hybrid_SOMP(:,i)= R_hybrid_SOMP(:,i)+SE_BER(2, H, G_SOMP, W_RF_SOMP, D2, F_RF_SOMP, channel_power, snr);
            R_hybrid_SOMP1(:,i) = R_hybrid_SOMP1(:,i)+SE_BER(2, H_temp, G_SOMP1, W_RF_SOMP_1, D21, F_RF_SOMP_1, channel_power, snr);
            R_hybrid_SOMP2(:,i)=  R_hybrid_SOMP2(:,i)+SE_BER(2, H, G_SOMP2, W_RF_SOMP_2, D22, F_RF_SOMP_2, channel_power, snr);
            R_hybrid_SOMP3 (:,i)=  R_hybrid_SOMP3 (:,i)+SE_BER(2, H_temp, G_SOMP3, W_RF_SOMP_3, D23, F_RF_SOMP_3, channel_power, snr);
            R_hybrid_SOMP4 (:,i)=  R_hybrid_SOMP4 (:,i)+SE_BER(2, H, G_SOMP4, W_RF_SOMP_4, D24, F_RF_SOMP_4, channel_power, snr);
            R_hybrid_SOMP5 (:,i)=  R_hybrid_SOMP5 (:,i)+SE_BER(2, H_temp, G_SOMP5, W_RF_SOMP_5, D25, F_RF_SOMP_5, channel_power, snr);
            R_hybrid_SOMP6 (:,i)=  R_hybrid_SOMP6 (:,i)+SE_BER(2, H, G_SOMP6, W_RF_SOMP_6, D26, F_RF_SOMP_6, channel_power, snr);
            R_hybrid_SOMP7(:,i)=  R_hybrid_SOMP7 (:,i)+SE_BER(2, H_temp, G_SOMP7, W_RF_SOMP_7, D27, F_RF_SOMP_7, channel_power, snr);
       fprintf('SNR=%d dB,iter=%d\n',SNR_dB(i),i_iter );
    end
      toc

end

R_hybrid_PCA=R_hybrid_PCA/Nsym;
R_hybrid_SOMP=R_hybrid_SOMP/Nsym;
R_hybrid_SOMP1=R_hybrid_SOMP1/Nsym;
R_hybrid_SOMP2=R_hybrid_SOMP2/Nsym;
R_hybrid_SOMP3=R_hybrid_SOMP3/Nsym;
R_hybrid_SOMP4=R_hybrid_SOMP4/Nsym;
R_hybrid_SOMP5=R_hybrid_SOMP5/Nsym;
R_hybrid_SOMP6=R_hybrid_SOMP6/Nsym;
R_hybrid_SOMP7=R_hybrid_SOMP7/Nsym;
            
%%  BER PLOT
% plot(SNR_dB,smooth(R_digital_GMD),'k-.','Linewidth',1.5)
% hold on
plot(SNR_dB, R_hybrid_PCA ,'b-o','Linewidth',1.5)
hold on
plot(SNR_dB,R_hybrid_SOMP,'g-s','Linewidth',1.5)
hold on
plot(SNR_dB,R_hybrid_SOMP1,'g-v','Linewidth',1.5)
hold on
plot(SNR_dB, R_hybrid_SOMP2,'r-+','Linewidth',1.5)
hold on
plot(SNR_dB,R_hybrid_SOMP3,'r-d','Linewidth',1.5)
hold on
plot(SNR_dB, R_hybrid_SOMP4,'k->','Linewidth',1.5)
hold on
plot(SNR_dB,R_hybrid_SOMP5,'k-^','Linewidth',1.5)
hold on
plot(SNR_dB, R_hybrid_SOMP6,'y-*','Linewidth',1.5)
hold on
plot(SNR_dB,R_hybrid_SOMP7,'y-p','Linewidth',1.5)
hold on
xlabel('SNR (dB)')
ylabel('Spectrum Efficiency(bps/Hz)')
legend('PCA with RIS','Ideal SOMP with RIS','Ideal SOMP without RIS','\rho=1 SOMP with RIS ','\rho=1 SOMP without RIS','\rho=2 SOMP with RIS ','\rho=2 SOMP without RIS','\rho=3 SOMP with RIS ','\rho=3 SOMP without RIS');
grid on


