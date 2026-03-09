function [SAT, UAV, UE] = FUNC_Create_Model(Area, numSAT, numUAV, numUE, zs, zu)
x_area = Area(1); y_area = Area(2);
numUAV_s = floor(numUAV/numSAT);
numUE_s = floor(numUE/numSAT);
% Create SATs
% Firstly, assume numSAT = 4
SAT = zeros(3, 4);
SAT(1,:) = [x_area/4 x_area*3/4 x_area*3/4 x_area/4];
SAT(2,:) = [y_area/4 y_area/4 y_area*3/4 y_area*3/4];
SAT(3,:) = [zs zs zs zs];

% Create UAVs
UAV = [];
UAV_1 = zu*ones(3, numUAV_s);
UAV_1(1,:) = x_area/2*rand(1,numUAV_s);
UAV_1(2,:) = y_area/2*rand(1,numUAV_s);
UAV = [UAV UAV_1];
UAV_2 = zu*ones(3, numUAV_s);
UAV_2(1,:) = x_area/2 + x_area/2*rand(1,numUAV_s);
UAV_2(2,:) = y_area/2*rand(1,numUAV_s);
UAV = [UAV UAV_2];
UAV_3 = zu*ones(3, numUAV_s);
UAV_3(1,:) = x_area/2 + x_area/2*rand(1,numUAV_s);
UAV_3(2,:) = y_area/2 + y_area/2*rand(1,numUAV_s);
UAV = [UAV UAV_3];
UAV_4 = zu*ones(3, numUAV_s);
UAV_4(1,:) = x_area/2*rand(1,numUAV_s);
UAV_4(2,:) = y_area/2 + y_area/2*rand(1,numUAV_s);
UAV = [UAV UAV_4];

% Create UEs
% UE = zeros(3, numUE);
% UE(1,:) = x_area*rand(1,numUE);
% UE(2,:) = y_area*rand(1,numUE);
UE = [];
UE_1 = zeros(3, numUE_s);
UE_1(1,:) = x_area/2*rand(1,numUE_s);
UE_1(2,:) = y_area/2*rand(1,numUE_s);
UE = [UE UE_1];
UE_2 = zeros(3, numUE_s);
UE_2(1,:) = x_area/2 + x_area/2*rand(1,numUE_s);
UE_2(2,:) = y_area/2*rand(1,numUE_s);
UE = [UE UE_2];
UE_3 = zeros(3, numUE_s);
UE_3(1,:) = x_area/2 + x_area/2*rand(1,numUE_s);
UE_3(2,:) = y_area/2 + y_area/2*rand(1,numUE_s);
UE = [UE UE_3];
UE_4 = zeros(3, numUE_s);
UE_4(1,:) = x_area/2*rand(1,numUE_s);
UE_4(2,:) = y_area/2 + y_area/2*rand(1,numUE_s);
UE = [UE UE_4];

end