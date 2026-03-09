clc
clear
close all

% Scenario
N = 100; % The number of radiation feeds/antennas
B = 25; % The maximum beams per SAT
numSAT = 4; % The number of SATs S = 4
numUAVeachSAT = 4;
numUAV = numUAVeachSAT*numSAT; % The number of UAVs U = 16 = 4 * 4
numUEeachUAV = 10; % average
numUE = numUAV*numUEeachUAV; % The number of UEs K = 80
numF = 30; % The number of possible packages F = 30
numM = 5; % The cache capacity M = 5
N_U = 20; % The maximum number os UEs each UAV can serve
v_light = 3e8; % The velocity of light
Q = 1e3; % The size of one package
% shadow_deviation = 4; % standard deviation shadowing 8dB, NLOS

PmaxSAT = 10^(50/10 - 3); % The maximum power of each SAT 50 dBm
PmaxUAV = 5; % 5W %10^(20/10 - 3); % The maximum power of each UAV 20 dBm
% P0 = 10^(10/10); % Pcircuit = 10 dBW
fc = 2e9; % 2GHz
BW = 20*1e6; % 10MHz
noiseVariancedBm = -174 + 10*log10(BW);
Noise_var=db2pow(noiseVariancedBm-30);
%Noise_var = 1;
Area = [10e3 10e3]; % 20kmx20km
zs = 780e3; % The altitude of SATs
zu = 0.5e3; % The altitude of UAVs


% satcom = zeros(3,100);
% for i = 1:100
%% Create model
MODEL = struct();
[SAT, UAV, UE] = FUNC_Create_Model(Area, numSAT, numUAV, numUE, zs, zu);
MODEL.SAT = SAT; MODEL.UAV = UAV; MODEL.UE = UE;
% FUNC_Plot_model_2D(MODEL, Area, "normal");
MODEL.H_SU = FUNC_Path_SAT_UAV(MODEL, N);
MODEL.H_UU = FUNC_Path_UAV_UAV(MODEL);
MODEL.H_UK = FUNC_Path_UAV_UE(MODEL, fc);

% Precoding
H_norm2 = sqrt(abs(sum( conj(MODEL.H_SU).*MODEL.H_SU, 3)));
MODEL.PrecodingS = conj(MODEL.H_SU)./ repmat(H_norm2,1,1,N);


%% Initialize RP (requested packages), A, B, and P
[RP, A, B, P_sat, P_uav_uav, P_uav_ue] = FUNC_INITIALIZE(numSAT, numUAV, numUE, N_U, ...
    numF, numM, PmaxSAT, PmaxUAV);
MODEL.RP = RP;
% save('model.mat');
% k = 1;
% f = RP(k);
% [tk, case_tk] = FUNCLatency(MODEL, k, f, A, B, P_sat, P_uav_uav, P_uav_ue,...
%               Area, Q, Noise_var, BW);
% satcom(case_tk,i) = tk;  
% end


%% Begin processing
A_random = FUNC_RandomClustering(numSAT, numUAV, numUE, N_U);
A_nearest = FUNC_NearestClustering(MODEL, numSAT, numUAV, numUE, N_U);
B_random = FUNC_RandomCache(numUAV, numM, numF);
[t_total, T] = func_ComputeSumLatency(MODEL, A_random, B_random,...
              Area, Q, Noise_var, BW, numUAVeachSAT, numUEeachUAV, numSAT, numUAV, numUE, PmaxSAT, PmaxUAV);

% tic
% A = FUNC_GameTheory_Clustering(MODEL, A, B, numSAT, numUAV,...
%     numUE, numUAVeachSAT, numUEeachUAV, N_U, PmaxSAT, PmaxUAV, Area, Q, Noise_var, BW);
% toc
% tic
% P = 40; % #elements in the population
% Pc = 0.7; % The probability of crossover
% Pm = 0.3; % The probability of mutation
% s = 1;
% [B_full, BestFit] = FUNC_GA(P, Pc, Pm,  s, MODEL, A, B, Area, Q, Noise_var,...
%     BW, numUAVeachSAT, numUEeachUAV, numM, numF, numSAT, numUAV, numUE, PmaxSAT, PmaxUAV);
% toc


%% Save to pdf
function save_to_pdf(num)
%Position plot at left hand corner with width 5 and height 5.
set(figure(num), 'PaperPosition', [0 0 5 5]);
%Set the paper to have width 5 and height 5.
set(figure(num), 'PaperSize', [5 5]); 
saveas(figure(num), 'GT_Clustering', 'pdf')
end
