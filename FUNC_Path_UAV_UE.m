function H = FUNC_Path_UAV_UE(MODEL, fc)
numUAV = size(MODEL.UAV,2); numUE = size(MODEL.UE,2);
% Distance from SATs to UAVs, row-SAT, column-UE
Distance = sqrt( (MODEL.UAV(1,:)'-MODEL.UE(1,:)).^2 + ...
                 (MODEL.UAV(2,:)'-MODEL.UE(2,:)).^2 + ...
                 (MODEL.UAV(3,:)'-MODEL.UE(3,:)).^2 );
             
H = zeros(numUAV, numUE);
for u = 1:numUAV
    for k = 1:numUE
        dist = Distance(u,k);
        z_u = MODEL.UAV(3,u);
        dist_xy = sqrt(dist^2 - z_u^2);
        H(u,k) = FUNC_Path_1UAV_1UE(1, 1, dist, dist_xy, z_u, fc);
    end
end
end


function Channel_1BS_1UE = FUNC_Path_1UAV_1UE(T, R, dist, dist_xy, z_u, fc)
pl_exponent = 2;   % path-loss exponent
% fc = 20*10^9;
c = 3*10^8;
a = 9.61; b = 0.16; % ATG channel parameters
eta_los = 1;
eta_nlos = 20;

PL_mk = 10*pl_exponent*log10(4*pi*fc*dist/c);
Pr_los = 1/(1+a*exp(-b*(atan(z_u/dist_xy) - a)));
Pr_nlos = 1 - Pr_los;
patlossSha_mac = PL_mk + eta_los*Pr_los + eta_nlos*Pr_nlos;

% Small-scale fading
Small_scale = (randn(R,T)+1i*randn(R,T))/sqrt(2);

% Combine large-scale fading and small-scale fading
Channel_1BS_1UE = sqrt(1/(10^(patlossSha_mac/10))).*Small_scale;
end
