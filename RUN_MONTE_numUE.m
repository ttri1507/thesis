clc
clear
close all

% Scenario
N = 100; % The number of radiation feeds/antennas
B = 25; % The maximum beams per SAT
numSAT = 4; % The number of SATs M = 4

Pmax = 10^(20/10); % The maximum power of each SAT 20 dBW
P0 = 10^(10/10); % Pcircuit = 10 dBW
Bw = 500*1e6; % 500MHz
% noiseVariancedBm = -174 + 10*log10(Bw) + 9;
% Noise_var=db2pow(noiseVariancedBm-30);
Noise_var = 1;
f = 20e9; % 20GHz
Area = [2000e3 2000e3]; % 2000kmx2000km
zm = 1000e3; % The altitude of SATs
t_thres = 2000e3/3e8; % The distance must be less than 1546km
Rmin = 0.2; % xxxx the minimum data rate each UE

% loop
numMonte = 100;

NumUE = [50 60 70 80 90 100];

EE_total = zeros(4,length(NumUE));

for numScen = 1:length(NumUE)
fprintf('SCENARIO NUMBER %d \n', numScen);
numCurrent = 1;
EE_Monte = zeros(4, numMonte);
% Command this
while(numCurrent <= numMonte)
    
% Create model
MODEL = struct();
[SAT, UE] = FUNC_Create_SAT_UE(Area, numSAT, NumUE(numScen), zm);
MODEL.SAT = SAT; MODEL.UE = UE;
%FUNC_Plot_model_2D(MODEL, Area, "normal");
MODEL.Channel = FUNC_PATH_SAT(MODEL, N, Bw, f);

% Clustering based Coalition game
theta = 0.5; % Trade-off between latency and channel gain
% start_RC = tic;
[MODEL_RC] = FUNC_randomClustering(MODEL,B,t_thres,theta,Area);
% time_RC = toc(start_RC);
% figure(1),FUNC_Plot_model_2D(MODEL_RC, Area, "cluster")

IdClusUE_RC = MODEL_RC.IdClusUE; NumUE_in_SAT_RC = MODEL_RC.NumUE_in_SAT;
start_GTC = tic;
[MODEL_GTC] = FUNC_CoalitionGame(MODEL,B,t_thres,theta,Area,IdClusUE_RC,NumUE_in_SAT_RC);
time_GTC = toc(start_GTC);
% figure(2),FUNC_Plot_model_2D(MODEL_GTC, Area, "cluster")

% Precoding
H_norm2 = sqrt(abs(sum( conj(MODEL.Channel).*MODEL.Channel, 3)));
% precoding for all UEs
Precoding = conj(MODEL.Channel)./ repmat(H_norm2,1,1,N);

% Power allocation
% Rho square
Rho = abs(sum(MODEL.Channel.*Precoding,3)).^2;
sigma_u = Noise_var; % sigma square
Pow_overline = sigma_u*(2^Rmin - 1)./Rho;
Beta = 1 + Pow_overline.*Rho/sigma_u;

% %%%%%%%%%%%% GTC-PA
NumUE_in_SAT = MODEL_GTC.NumUE_in_SAT;
IdClusUE = MODEL_GTC.IdClusUE;
lambda_ref = 1; tau_ref2 = 1;
tau_underline2 = 1e-10; tau_overline2 = 1000;
SumP_overline = zeros(1, numSAT);
for m = 1:numSAT
    Um = find(IdClusUE==m);
    SumP_overline(m) = sum(Pow_overline(m,Um));
end
Pmax_overline = Pmax - SumP_overline;
Pow_cir_overline = P0 + SumP_overline;
[sumrate_GTCPA, sumpow_GTCPA, obj_GTCPA, Pow_GTCPA, Rate_GTCPA] = FUNC_multiSAT_OPT(NumUE_in_SAT, IdClusUE, lambda_ref, Beta, ...
            Rho, sigma_u, Pmax_overline, Pow_cir_overline, tau_ref2, tau_underline2, tau_overline2);
EE_GTCPA = sumrate_GTCPA/(sum(Pow_cir_overline)+sumpow_GTCPA);

% %%%%%%%%%%%% RC-PA
NumUE_in_SAT = MODEL_RC.NumUE_in_SAT;
IdClusUE = MODEL_RC.IdClusUE;
lambda_ref = 1; tau_ref2 = 1;
tau_underline2 = 1e-10; tau_overline2 = 1000;
SumP_overline = zeros(1, numSAT);
for m = 1:numSAT
    Um = find(IdClusUE==m);
    SumP_overline(m) = sum(Pow_overline(m,Um));
end
Pmax_overline = Pmax - SumP_overline;
Pow_cir_overline = P0 + SumP_overline;
[sumrate_RCPA, sumpow_RCPA, obj_RCPA, Pow_RCPA, Rate_RCPA] = FUNC_multiSAT_OPT(NumUE_in_SAT, IdClusUE, lambda_ref, Beta, ...
            Rho, sigma_u, Pmax_overline, Pow_cir_overline, tau_ref2, tau_underline2, tau_overline2);
EE_RCPA = sumrate_RCPA/(sum(Pow_cir_overline)+sumpow_RCPA);
        
% %%%%%%%%%%%% GTC-EP
NumUE_in_SAT = MODEL_GTC.NumUE_in_SAT; IdClusUE = MODEL_GTC.IdClusUE;
P_eachUE = Pmax./NumUE_in_SAT;
Pow_GTCEP = zeros(numSAT,NumUE(numScen));
for u = 1:NumUE(numScen)
    m = IdClusUE(u);
    Pow_GTCEP(m,u) = P_eachUE(m);
end
[sumrate_GTCEP, sumpow_GTCEP, Rate_GTCEP] = FUNC_Compute_values(NumUE_in_SAT, IdClusUE, Rho, sigma_u, Pow_GTCEP);
EE_GTCEP = sumrate_GTCEP/(sumpow_GTCEP + numSAT*P0);

% %%%%%%%%%%%% RC-EP
NumUE_in_SAT = MODEL_RC.NumUE_in_SAT; IdClusUE = MODEL_RC.IdClusUE;
P_eachUE = Pmax./NumUE_in_SAT;
Pow_RCEP = zeros(numSAT,NumUE(numScen));
for u = 1:NumUE(numScen)
    m = IdClusUE(u);
    Pow_RCEP(m,u) = P_eachUE(m);
end
[sumrate_RCEP, sumpow_RCEP, Rate_RCEP] = FUNC_Compute_values(NumUE_in_SAT, IdClusUE, Rho, sigma_u, Pow_RCEP);
EE_RCEP = sumrate_RCEP/(sumpow_RCEP + numSAT*P0);

EE_Monte(1,numCurrent) = EE_GTCPA; EE_Monte(2,numCurrent) = EE_RCPA;
EE_Monte(3,numCurrent) = EE_GTCEP; EE_Monte(4,numCurrent) = EE_RCEP;
numCurrent = numCurrent + 1;
end %end while
EE_total(:,numScen) = sum(EE_Monte,2)/numMonte;
end

% Show the result
X = NumUE;
semilogy(X,EE_total(1,:),'r--^','markersize',4,'Linewidth',1);
hold on; grid on;
semilogy(X,EE_total(2,:),'b--^','markersize',4,'Linewidth',1);
semilogy(X,EE_total(3,:),'r-*','markersize',4,'Linewidth',1);
semilogy(X,EE_total(4,:),'b-*','markersize',4,'Linewidth',1);
legend('EE_GTCPA','EE_RCPA','EE_GTCEP','EE_RCEP')
xlabel('Number of UEs')
ylabel('EE (bits/Joule/Hz)')

%% Save EE_total to txt file
% Data = [X', EE_total'];
% fileID = fopen('EE-UE.txt','w');
% fprintf(fileID,'%20s %20s %20s %20s %20s\n','X','EE_GTCPA','EE_RCPA','EE_GTCEP','EE_RCEP');
% for i = 1:length(NumUE)
%     fprintf(fileID,'%20.5f %20.5f %20.5f %20.5f %20.5f\n',Data(i,:));
% end
% fclose(fileID);
