clc
clear
close all

% Scenario
N = 100; % The number of radiation feeds/antennas
B = 25; % The maximum beams per SAT
numSAT = 4; % The number of SATs S = 4

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
numMonte = 20;

NumUAVeachSAT = [2 4 6 8 10];
NumUAV = NumUAVeachSAT*numSAT; % The number of UAVs U = 16 = 4 * 4
numUEeachUAV = 5;
NumUE = NumUAV*numUEeachUAV;

Time_total = zeros(2,length(NumUAVeachSAT));

for numScen = 1:length(NumUAVeachSAT)
fprintf('SCENARIO NUMBER %d \n', numScen);
numCurrent = 1;
Time_Monte = zeros(2, numMonte);
% Command this
while(numCurrent <= numMonte)
MODEL = struct();
[SAT, UAV, UE] = FUNC_Create_Model(Area, numSAT, NumUAV(numScen), NumUE(numScen), zs, zu);
MODEL.SAT = SAT; MODEL.UAV = UAV; MODEL.UE = UE;
% FUNC_Plot_model_2D(MODEL, Area, "normal");
MODEL.H_SU = FUNC_Path_SAT_UAV(MODEL, N);
MODEL.H_UU = FUNC_Path_UAV_UAV(MODEL);
MODEL.H_UK = FUNC_Path_UAV_UE(MODEL, fc);

% Precoding
H_norm2 = sqrt(abs(sum( conj(MODEL.H_SU).*MODEL.H_SU, 3)));
MODEL.PrecodingS = conj(MODEL.H_SU)./ repmat(H_norm2,1,1,N);


% Initialize RP (requested packages), A, B, and P
[RP, A_ini, B_ini, P_sat_ini, P_uav_uav_ini, P_uav_ue_ini] = FUNC_INITIALIZE(numSAT, NumUAV(numScen), NumUE(numScen), N_U, ...
    numF, numM, PmaxSAT, PmaxUAV);
MODEL.RP = RP;
% save('model.mat');

% COMPARE
P = 20; % #elements in the population
Pc = 0.8; % The probability of crossover
Pm = 0.2; % The probability of mutation

%%%%%%%%%%%%%%%%% 1. Game Theory and Genetic Algorithm (GTGA) %%%%%%%%%%%%%%%%%
A_GT = A_ini; B_GA = B_ini; t_1 = 1e6;
start_GT = tic;
A_GT = FUNC_GameTheory_Clustering(MODEL, A_GT, B_GA, numSAT, NumUAV(numScen),...
        NumUE(numScen), NumUAVeachSAT(numScen), numUEeachUAV, N_U, PmaxSAT, PmaxUAV, Area, Q, Noise_var, BW);
time_GT = toc(start_GT);

time_GA = 0;
for s = 1:4
    start_GA = tic;
    [B_GA, BestFit] = FUNC_GA(P, Pc, Pm,  s, MODEL, A_GT, B_GA, Area, Q, Noise_var,...
         BW, NumUAVeachSAT(numScen), numUEeachUAV, numM, numF, numSAT, NumUAV(numScen), NumUE(numScen), PmaxSAT, PmaxUAV);
    time_GA = time_GA + toc(start_GA);
end
time_GA = time_GA/4;
[t_current, T] = func_ComputeSumLatency(MODEL, A_GT, B_GA,...
                  Area, Q, Noise_var, BW, NumUAVeachSAT(numScen), numUEeachUAV, numSAT, NumUAV(numScen), NumUE(numScen), PmaxSAT, PmaxUAV);
%     fprintf('t_total = %f\n', t_1);

Time_Monte(1,numCurrent) = time_GT;
Time_Monte(2,numCurrent) = time_GA;

numCurrent = numCurrent + 1;
end %end while
Time_total(:,numScen) = sum(Time_Monte,2)/numMonte;
end

% Show the result
X = NumUE;
semilogy(X,Time_total(1,:),'r--^','markersize',4,'Linewidth',1);
hold on; grid on;
semilogy(X,Time_total(2,:),'b--^','markersize',4,'Linewidth',1);
legend('Time_GT','Time_GA')
xlabel('Number of UEs')
ylabel('Executive time (s)')

% % Save Time_total to txt file
% Data = [X', Time_total'];
% fileID = fopen('Time-UAV-UE.txt','w');
% fprintf(fileID,'%20s %20s %20s\n','X','Time_GT','Time_GA');
% for i = 1:length(NumUE)
%     fprintf(fileID,'%20.5f %20.5f %20.5f\n',Data(i,:));
% end
% fclose(fileID);