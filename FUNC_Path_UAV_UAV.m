function H = FUNC_Path_UAV_UAV(MODEL)
numUAV = size(MODEL.UAV,2);
% Distance from SATs to UAVs, row-SAT, column-UE
Distance = sqrt( (MODEL.UAV(1,:)'-MODEL.UAV(1,:)).^2 + ...
                 (MODEL.UAV(2,:)'-MODEL.UAV(2,:)).^2 + ...
                 (MODEL.UAV(3,:)'-MODEL.UAV(3,:)).^2 );
             
H = zeros(numUAV, numUAV);
for up = 1:numUAV
    for u = 1:numUAV
        dist = Distance(up,u);
        H(up,u) = FUNC_Path_1UAV_1UAV(1, 1, dist);
    end
end
end


function Channel_1BS_1UE = FUNC_Path_1UAV_1UAV(T, R, dist)
pl_freespace = 2;
beta_0 = 10^(-30/10); % Channel power gain at the reference distance  
UAV_UAV_path = beta_0*dist^(-pl_freespace);

% Small-scale fading
Small_scale = (randn(R,T)+1i*randn(R,T))/sqrt(2);

% Combine large-scale fading and small-scale fading
Channel_1BS_1UE = sqrt(UAV_UAV_path).*Small_scale;
end
