% H = SAT x UE x N
% N is the number of antennas per SAT
function H = FUNC_Path_SAT_UAV(MODEL, N)
numSAT = size(MODEL.SAT,2); numUAV = size(MODEL.UAV,2); numUE = size(MODEL.UE,2);
% Distance from SATs to UAVs, row-SAT, column-UE
Distance = sqrt( (MODEL.SAT(1,:)'-MODEL.UAV(1,:)).^2 + ...
                 (MODEL.SAT(2,:)'-MODEL.UAV(2,:)).^2 + ...
                 (MODEL.SAT(3,:)'-MODEL.UAV(3,:)).^2 );
             
H = zeros(numSAT, numUAV, N);
for m = 1:numSAT
    for u = 1:numUAV
        H(m,u,:) = func_path_1sat_1ss(N, 1, Distance(m,u));
    end
end
end


% N: the number of radiation elements
% d_mu: the distance from SAT to sensor u
function Channel_1SAT_1UAV = func_path_1sat_1ss(T, R, dist)
% omega = om (LoS) = d, deta = b (scatter) = sigma, epsilon = m (Nakagami)
% ref: https://github.com/Jonathan-Browning/Shadowed-Rician-Fading-Matlab/blob/main/src/Class/ShadowedRice.m
% https://github.com/Jonathan-Browning/Shadowed-Rician-Fading-Matlab/blob/main/src/Class/ShadowedRice.m
% Trajectory Design and Link Selection in UAV-assisted Hybrid Satellite-terrestrial Network
pl_freespace = 2;
beta_0 = 10^(-30/10); % Channel power gain at the reference distance  
SAT_UAV_path = beta_0*dist^(-pl_freespace);

% SR fading
om = 0.0005; b = 0.063; m = 2;
% Small_scale = randn(R,T)+1i*randn(R,T);
naka_m = sqrt(gamrnd(m,1/m));
Phi = 2*pi*rand(T,R);
p = sqrt(om)*cos(Phi); q = sqrt(om)*sin(Phi);
sigma = sqrt(b);
x = normrnd(naka_m*p, sigma); y = normrnd(naka_m*q, sigma);
sr_fading = x + 1i*y;

% Combine large-scale fading and small-scale fading
Channel_1SAT_1UAV = sqrt(SAT_UAV_path.*sr_fading);%.*Small_scale;
end
