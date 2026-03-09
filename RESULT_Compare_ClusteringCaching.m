clc
clear
close all

% Scenario
N = 100; % The number of radiation feeds/antennas
B = 25; % The maximum beams per SAT
numSAT = 4; % The number of SATs S = 4
numUAVeachSAT = 4;
numUAV = numUAVeachSAT*numSAT; % The number of UAVs U = 16 = 4 * 4

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

% MONTE CARLO
numMonte = 100;

NumUEeachUAV = [5 6 7 8 9 10]; %
NumUE = numUAV.*NumUEeachUAV; %

T_total = zeros(6,length(NumUEeachUAV));

for numScen = 1:length(NumUEeachUAV)
fprintf('SCENARIO NUMBER %d \n', numScen);
numCurrent = 1;
T_Monte = zeros(6, numMonte);
% Command this
while(numCurrent <= numMonte)
MODEL = struct();
[SAT, UAV, UE] = FUNC_Create_Model(Area, numSAT, numUAV, NumUE(numScen), zs, zu);
MODEL.SAT = SAT; MODEL.UAV = UAV; MODEL.UE = UE;
% FUNC_Plot_model_2D(MODEL, Area, "normal");
MODEL.H_SU = FUNC_Path_SAT_UAV(MODEL, N);
MODEL.H_UU = FUNC_Path_UAV_UAV(MODEL);
MODEL.H_UK = FUNC_Path_UAV_UE(MODEL, fc);

% Precoding
H_norm2 = sqrt(abs(sum( conj(MODEL.H_SU).*MODEL.H_SU, 3)));
MODEL.PrecodingS = conj(MODEL.H_SU)./ repmat(H_norm2,1,1,N);


% Initialize RP (requested packages), A, B, and P
[RP, A_ini, B_ini, P_sat_ini, P_uav_uav_ini, P_uav_ue_ini] = FUNC_INITIALIZE(numSAT, numUAV, NumUE(numScen), N_U, ...
    numF, numM, PmaxSAT, PmaxUAV);
MODEL.RP = RP;
% save('model.mat');

% COMPARE
P = 20; % #elements in the population
Pc = 0.8; % The probability of crossover
Pm = 0.2; % The probability of mutation

% %%%%%%%%%%%%%%%%% 1. Game Theory and Genetic Algorithm (GTGA) %%%%%%%%%%%%%%%%%
% A_GT = A_ini; B_GA = B_ini; t_1 = 1e6;
% for loop1 = 1:5
%     A_GT = FUNC_GameTheory_Clustering(MODEL, A_GT, B_GA, numSAT, numUAV,...
%         NumUE(numScen), numUAVeachSAT, NumUEeachUAV(numScen), N_U, PmaxSAT, PmaxUAV, Area, Q, Noise_var, BW);
%     for s = 1:4
%         [B_GA, BestFit] = FUNC_GA(P, Pc, Pm,  s, MODEL, A_GT, B_GA, Area, Q, Noise_var,...
%             BW, numUAVeachSAT, NumUEeachUAV(numScen), numM, numF, numSAT, numUAV, NumUE(numScen), PmaxSAT, PmaxUAV);
%     end
%     [t_1_current, T_1] = func_ComputeSumLatency(MODEL, A_GT, B_GA,...
%                   Area, Q, Noise_var, BW, numUAVeachSAT, NumUEeachUAV(numScen), numSAT, numUAV, NumUE(numScen), PmaxSAT, PmaxUAV);
% %     fprintf('t_total = %f\n', t_1);
%     if(t_1_current<t_1)
%         t_1 = t_1_current;
%     end
% end
% 
% %%%%%%%%%%%%%%%%% 2. Random Clustering and Genetic Algorithm (RCGA) %%%%%%%%%%%%%%%%%
% A_random = FUNC_RandomClustering(numSAT, numUAV, NumUE(numScen), N_U);
% B_GA = B_ini;
% for s = 1:4
%     [B_GA, BestFit] = FUNC_GA(P, Pc, Pm,  s, MODEL, A_random, B_GA, Area, Q, Noise_var,...
%         BW, numUAVeachSAT, NumUEeachUAV(numScen), numM, numF, numSAT, numUAV, NumUE(numScen), PmaxSAT, PmaxUAV);
% end
% [t_2, T_2] = func_ComputeSumLatency(MODEL, A_random, B_GA,...
%               Area, Q, Noise_var, BW, numUAVeachSAT, NumUEeachUAV(numScen), numSAT, numUAV, NumUE(numScen), PmaxSAT, PmaxUAV);
% %%%%%%%%%%%%%%%%% 3. Nearest Clustering and Genetic Algorithm (NCGA)%%%%%%%%%%%%%%%%%
% A_nearest = FUNC_NearestClustering(MODEL, numSAT, numUAV, NumUE(numScen), N_U);
% B_GA = B_ini;
% for s = 1:4
%     [B_GA, BestFit] = FUNC_GA(P, Pc, Pm,  s, MODEL, A_nearest, B_GA, Area, Q, Noise_var,...
%         BW, numUAVeachSAT, NumUEeachUAV(numScen), numM, numF, numSAT, numUAV, NumUE(numScen), PmaxSAT, PmaxUAV);
% end
% [t_3, T_3] = func_ComputeSumLatency(MODEL, A_nearest, B_GA,...
%               Area, Q, Noise_var, BW, numUAVeachSAT, NumUEeachUAV(numScen), numSAT, numUAV, NumUE(numScen), PmaxSAT, PmaxUAV);
%%%%%%%%%%%%%%%%% 4. Game Theory and Random Cache Placement (GTRP) %%%%%%%%%%%%%%%%%
A_GT = FUNC_GameTheory_Clustering(MODEL, A_ini, B_ini, numSAT, numUAV,...
    NumUE(numScen), numUAVeachSAT, NumUEeachUAV(numScen), N_U, PmaxSAT, PmaxUAV, Area, Q, Noise_var, BW);
B_random = FUNC_RandomCache(numUAV, numM, numF);
[t_4, T_4] = func_ComputeSumLatency(MODEL, A_GT, B_random,...
              Area, Q, Noise_var, BW, numUAVeachSAT, NumUEeachUAV(numScen), numSAT, numUAV, NumUE(numScen), PmaxSAT, PmaxUAV);
%%%%%%%%%%%%%%%%% 5. Random Clustering and Random Cache Placement (RCRP) %%%%%%%%%%%%%%%%%
A_random = FUNC_RandomClustering(numSAT, numUAV, NumUE(numScen), N_U);
B_random = FUNC_RandomCache(numUAV, numM, numF);
[t_5, T_5] = func_ComputeSumLatency(MODEL, A_random, B_random,...
              Area, Q, Noise_var, BW, numUAVeachSAT, NumUEeachUAV(numScen), numSAT, numUAV, NumUE(numScen), PmaxSAT, PmaxUAV);
%%%%%%%%%%%%%%%%% 6. Nearest Clustering and Random Cache Placement (NCRP) %%%%%%%%%%%%%%%%%
A_nearest = FUNC_NearestClustering(MODEL, numSAT, numUAV, NumUE(numScen), N_U);
B_random = FUNC_RandomCache(numUAV, numM, numF);
[t_6, T_6] = func_ComputeSumLatency(MODEL, A_nearest, B_random,...
              Area, Q, Noise_var, BW, numUAVeachSAT, NumUEeachUAV(numScen), numSAT, numUAV, NumUE(numScen), PmaxSAT, PmaxUAV);
t_1 = 0; t_2 = 0; t_3 = 0;
T_Monte(1,numCurrent) = t_1; T_Monte(2,numCurrent) = t_2;
T_Monte(3,numCurrent) = t_3; T_Monte(4,numCurrent) = t_4;
T_Monte(5,numCurrent) = t_5; T_Monte(6,numCurrent) = t_6;
numCurrent = numCurrent + 1;
end %end while
T_total(:,numScen) = sum(T_Monte,2)/numMonte;
end

% Show the result
X = NumUE;
semilogy(X,T_total(1,:),'r--^','markersize',4,'Linewidth',1);
hold on; grid on;
semilogy(X,T_total(2,:),'r--o','markersize',4,'Linewidth',1);
semilogy(X,T_total(3,:),'r--*','markersize',4,'Linewidth',1);
semilogy(X,T_total(4,:),'b-^','markersize',4,'Linewidth',1);
semilogy(X,T_total(5,:),'b-o','markersize',4,'Linewidth',1);
semilogy(X,T_total(6,:),'b-*','markersize',4,'Linewidth',1);
legend('GTGA', 'RCGA', 'NCGA', 'GTRP', 'RCRP', 'NCRP')
xlabel('Number of GUs')
ylabel('Latency (s)')

%%Save T_total to txt file
% Data = [X', T_total'];
% fileID = fopen('CompareLatency.txt','w');
% fprintf(fileID,'%20s %20s %20s %20s %20s %20s %20s\n','X','GTGA', 'RCGA', 'NCGA', 'GTRP', 'RCRP', 'NCRP');
% for i = 1:length(NumUE)
%     fprintf(fileID,'%20.5f %20.5f %20.5f %20.5f %20.5f %20.5f %20.5f\n',Data(i,:));
% end
% fclose(fileID);