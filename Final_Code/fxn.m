function [rew] = fxn(p1)
% clc;
% clear all;
% close all;
%%
% p1=ones(10,1)*1;
p=zeros(10,1);
for i = 1:10
    p(i,1)= p(i,1)+p1(i);
end
selectSimulationSetup = 1;
nbrOfSetups = 1;

%Number of channel realizations per setup
nbrOfRealizations = 1000;

if selectSimulationSetup == 1
    
    %Number of APs per setup
    L = 50;
    
    %Number of antennas per AP
    N = 2;
    
elseif selectSimulationSetup == 2
    
    %Number of APs per setup
    L = 25;
    
    %Number of antennas per AP
    N = 4;
    
end

%Number of UEs in the network
K = size(p,1);

%Length of coherence block
tau_c = 200;

%Length of pilot sequences
tau_p = 10;


%Prepare to save simulation results
SE_scalable_MR_tot = zeros(K,nbrOfSetups);
SE_scalable_LP_MMSE_tot = zeros(K,nbrOfSetups);
SE_scalable_P_MMSE_tot = zeros(K,nbrOfSetups);
SE_scalable_MMSE_tot = zeros(K,nbrOfSetups);
SE_all_MR_tot = zeros(K,nbrOfSetups);
SE_all_LP_MMSE_tot = zeros(K,nbrOfSetups);
SE_all_P_MMSE_tot = zeros(K,nbrOfSetups);
SE_all_MMSE_tot = zeros(K,nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs and APs at random locations
    
    [gainOverNoisedB,R,pilotIndex,D] = generateSetup(L,K,N,tau_p,1);

    %Generate channel realizations with estimates and estimation
    %error correlation matrices
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    
    %% Proposed Scalable Cell-Free Massive MIMO
    
    %Compute SE using Propositions 1 and 2
    [SE_MR,SE_LP_MMSE,SE_P_MMSE,SE_MMSE] = functionComputeSE_uplink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex);
    
    %Save SE values
    SE_scalable_MR_tot(:,n) = SE_MR;
    SE_scalable_LP_MMSE_tot(:,n) = SE_LP_MMSE;
    SE_scalable_P_MMSE_tot(:,n) = SE_P_MMSE;
    SE_scalable_MMSE_tot(:,n) = SE_MMSE;
    
    
    %% Original Cell-Free Massive MIMO
    
    %Define the case when all APs serve all UEs
    D_all = ones(L,K);
    
    %Compute SE using Propositions 1 and 2
    [SE_MR,SE_LP_MMSE,SE_P_MMSE,SE_MMSE] = functionComputeSE_uplink(Hhat,H,D_all,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex);
    
    %Save SE values
%    SE_all_MR_tot(:,n) = SE_MR;
    SE_all_LP_MMSE_tot(:,n) = SE_LP_MMSE;
    SE_all_P_MMSE_tot(:,n) = SE_P_MMSE;
    SE_all_MMSE_tot(:,n) = SE_MMSE;
    
    
    %Remove large matrices at the end of analyzing this setup
    %clear Hhat H B C R;
    
end


%% Reward
%rew = sum(SE_scalable_LP_MMSE_tot(:))/sqrt(sum(p(:)));
rew = sum(SE_scalable_P_MMSE_tot(:))

%%
% figure;
% hold on;
% plot(M);
% xlabel('Iteration');
% ylabel('Reward');
% axis tight;
% title('Performance of RL model with time');

% plot(sort(SE_scalable_P_MMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
% xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'P-MMSE (Scalable)','LP-MMSE (Scalable)','MR (Scalable)'},'Interpreter','Latex','Location','SouthEast');
% xlim([0 10]);
