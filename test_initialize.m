clc
clear
load('model.mat')

numUAVeachSAT = numUAV/numSAT;
numUEeachSAT = numUE/numSAT;

% Initialize the requested packages from GUs
RP = randi(numF, [1,numUE]);

% Initialize Beta
B = zeros(numUAV,numF);
b = ones(numUAV,numM);
B(:,1:numM) = b;

% Initialize A: using diagonal line assign 1
A = zeros(numUAV,numUE);
k = 0; numLoop = 0; conBreak = false;
for s = 1:4
    Us = numUAVeachSAT*(s-1) + 1:numUAVeachSAT*s;
    Ks = numUEeachSAT*(s-1) + 1:numUEeachSAT*s;
    for k = 1:numUEeachSAT
        u = mod(k,numUAVeachSAT);
        if(u==0)
            u = numUAVeachSAT;
        end
        indUAV = Us(u);
        indUE = Ks(k);
        A(indUAV,indUE) = 1;
    end
end

% Initialize Power
P_sat = zeros(numSAT, numUAV);
p_eachSat = 0.5*PmaxSAT/numUAVeachSAT*ones(1,numUAVeachSAT);
for s = 1:4
    indUAVbegin = numUAVeachSAT*(s - 1) + 1;
    P_sat(s, indUAVbegin:indUAVbegin+numUAVeachSAT-1) = p_eachSat;
end
P_uav_ue = A*PmaxUAV*0.75/N_U; % 75% serve GUs, 25% for inter-UAV connections
P_uav_uav = zeros(numUAV,numUAV);
P_sub_uu = PmaxUAV*0.25/(numUAVeachSAT - 1)*ones(numUAVeachSAT,numUAVeachSAT);
for u = 1:numUAVeachSAT
    P_sub_uu(u,u) = 0;
end
for s = 1:4
    indUAVbegin = numUAVeachSAT*(s - 1) + 1;
    P_uav_uav(indUAVbegin:indUAVbegin+numUAVeachSAT-1,...
        indUAVbegin:indUAVbegin+numUAVeachSAT-1) = P_sub_uu;
end



